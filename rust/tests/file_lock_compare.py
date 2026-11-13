#!/usr/bin/env python3
"""Compare native installer locks with Python flock across real processes.

python file_lock_compare.py compare --probe BIN --work-root NEW_DIRECTORY
No production cache, lock or environment is opened. All test-owned processes
are reaped, including a deliberately killed holder; no unrelated PID is used.
"""
import argparse
import errno
import fcntl
import hashlib
import json
import os
from pathlib import Path
import selectors
import subprocess
import sys


def open_lock(path, kind):
    return path.open('a') if kind == 'append' else os.open(path, os.O_RDONLY)


def close_lock(handle):
    handle.close() if hasattr(handle, 'close') else os.close(handle)


def python_hold(args):
    handle = open_lock(args.path, args.open)
    try:
        print('opened', flush=True)
        flags = fcntl.LOCK_EX | (fcntl.LOCK_NB if args.wait == 'nonblocking' else 0)
        fcntl.flock(handle, flags)
        print('acquired', flush=True)
        action = sys.stdin.readline().rstrip('\r\n')
        if action == 'release':
            fcntl.flock(handle, fcntl.LOCK_UN)
        elif action == 'exit':
            os._exit(0)
        elif action not in ('drop', ''):
            raise ValueError('invalid lifetime command')
    finally:
        close_lock(handle)


def compare(args):
    root = args.work_root.absolute()
    root.mkdir(parents=True, exist_ok=False)
    processes = []
    cases = []
    failures = []

    def start(implementation, path, kind='append', wait='nonblocking'):
        command = ([str(args.probe)] if implementation == 'native' else
                   [sys.executable, str(Path(__file__).absolute()), 'hold'])
        command += ['--path', str(path), '--open', kind, '--wait', wait]
        child = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE, bufsize=0)
        processes.append(child)
        return child

    def line(child, timeout=10):
        # Read unbuffered bytes so a second already-buffered line cannot be
        # mistaken for an empty OS pipe by select().
        with selectors.DefaultSelector() as selector:
            selector.register(child.stdout, selectors.EVENT_READ)
            if not selector.select(timeout):
                return None
        return child.stdout.readline().decode().rstrip('\r\n')

    def holding(implementation, path, kind='append', wait='nonblocking'):
        child = start(implementation, path, kind, wait)
        assert line(child) == 'opened'
        assert line(child) == 'acquired'
        return child

    def finish(child, action='release'):
        child.stdin.write((action+'\n').encode())
        child.stdin.close()
        code = child.wait(timeout=10)
        assert code == 0, (code, child.stderr.read().decode())

    def blocked(path, kind):
        child = start('native', path, kind)
        output, errors = child.communicate(timeout=10)
        assert child.returncode == 1 and output == b'opened\n', (output, errors)
        # Rust's public TryLockError::WouldBlock deliberately carries no raw
        # errno. Its io::Error conversion preserves the failure category, not
        # Python's numeric errno/diagnostic spelling. Compare that category and
        # actual contention; never count a successful acquisition as equivalent.
        assert errors == b'errno=None; kind=WouldBlock\n', errors
        handle = open_lock(path, kind)
        try:
            try:
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                assert error.errno == errno.EWOULDBLOCK
            else:
                raise AssertionError('Python unexpectedly acquired held lock')
        finally:
            close_lock(handle)

    def record(name, procedure):
        first_child = len(processes)
        try:
            procedure()
        except Exception as error:
            failures.append({'case': name, 'error': repr(error)})
        else:
            cases.append(name)
        finally:
            # One failed case must not leave a holder contaminating aliases in
            # another case. Only children created by this case are closed.
            for child in processes[first_child:]:
                if child.poll() is None:
                    child.kill()
                child.wait(timeout=10)

    try:
        for kind in ('append', 'directory'):
            for holder in ('native', 'python'):
                for action in ('release', 'drop', 'exit', 'kill'):
                    name = f'{kind}-{holder}-{action}'
                    path = root / name
                    path.mkdir() if kind == 'directory' else path.write_bytes(b'preserved\n')
                    before = path.stat()

                    def lifecycle(path=path, kind=kind, holder=holder, action=action, before=before):
                        child = holding(holder, path, kind)
                        blocked(path, kind)
                        if action == 'kill':
                            child.kill()  # Only this test's live Popen child.
                            child.wait(timeout=10)
                        else:
                            finish(child, action)
                        finish(holding('native', path, kind))
                        finish(holding('python', path, kind))
                        after = path.stat()
                        assert (before.st_dev, before.st_ino) == (after.st_dev, after.st_ino)
                        if kind == 'append':
                            assert path.read_bytes() == b'preserved\n'
                        else:
                            assert not list(path.iterdir())
                    record(name, lifecycle)

            for holder, waiter in (('python', 'native'), ('native', 'python')):
                name = f'{kind}-{holder}-blocks-{waiter}'
                path = root / name
                path.mkdir() if kind == 'directory' else path.write_bytes(b'')

                def contention(path=path, kind=kind, holder=holder, waiter=waiter):
                    first = holding(holder, path, kind)
                    second = start(waiter, path, kind, 'blocking')
                    assert line(second) == 'opened'
                    assert line(second, 0.2) is None and second.poll() is None
                    finish(first)
                    assert line(second) == 'acquired'
                    finish(second)
                record(name, contention)

        path = root/'space and Unicode λ.lock'
        path.write_bytes(b'unchanged')
        for alias_kind in ('symlink', 'hardlink'):
            alias = root/(alias_kind+'.lock')
            alias.symlink_to(path) if alias_kind == 'symlink' else os.link(path, alias)

            def aliases(alias=alias):
                first = holding('native', path)
                blocked(alias, 'append')
                finish(first)
                finish(holding('native', alias))
                assert path.read_bytes() == b'unchanged'
            record('same-inode-'+alias_kind, aliases)

        def creation():
            original = os.umask(0o027)
            try:
                native = root/'native-created.lock'
                python = root/'python-created.lock'
                finish(holding('native', native))
                finish(holding('python', python))
                assert native.stat().st_mode & 0o777 == python.stat().st_mode & 0o777 == 0o640
                assert native.read_bytes() == python.read_bytes() == b''
            finally:
                os.umask(original)
        record('creation-umask-and-empty-file', creation)

        for name, path, kind in [('missing-parent', root/'absent'/'x', 'append'),
                                 ('missing-directory', root/'no-directory', 'directory'),
                                 ('append-directory', root, 'append')]:
            def opening_failure(path=path, kind=kind):
                try:
                    handle = open_lock(path, kind)
                except OSError as expected:
                    expected_errno = expected.errno
                else:
                    close_lock(handle)
                    raise AssertionError('Python unexpectedly opened negative fixture')
                child = start('native', path, kind)
                output, errors = child.communicate(timeout=10)
                assert child.returncode == 1 and not output
                assert f'errno=Some({expected_errno})' in errors.decode(), errors
            record(name, opening_failure)
    finally:
        for child in processes:
            if child.poll() is None:
                child.kill()
            child.wait(timeout=10)
            for pipe in (child.stdin, child.stdout, child.stderr):
                if pipe and not pipe.closed:
                    pipe.close()
    report = {'passed': not failures, 'cases_passed': len(cases), 'cases': cases,
              'mismatches': failures, 'processes_reaped': len(processes),
              'python': sys.version, 'probe_sha256': hashlib.sha256(args.probe.read_bytes()).hexdigest(),
              'scope': 'Same-host Unix flock interoperability, lifetime, opening errors, aliases and modes; no caller switch or cross-node filesystem claim',
              'failure_contract': 'Exact opening errno and lock outcome/category; Rust TryLockError normalizes contention and does not retain Python raw errno or diagnostic text'}
    (root/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'cases'}))
    if failures:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    p = commands.add_parser('compare')
    p.add_argument('--probe', type=Path, required=True)
    p.add_argument('--work-root', type=Path, required=True)
    p.set_defaults(function=compare)
    p = commands.add_parser('hold')
    p.add_argument('--path', type=Path, required=True)
    p.add_argument('--open', choices=('append', 'directory'), required=True)
    p.add_argument('--wait', choices=('blocking', 'nonblocking'), required=True)
    p.set_defaults(function=python_hold)
    args = parser.parse_args()
    args.function(args)


if __name__ == '__main__':
    main()
