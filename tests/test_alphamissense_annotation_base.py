"""Tiny real BGZF/index fixtures; VEP/conversion/CV are never scientifically run."""
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("alpha_base_test", ROOT / "scripts/alphamissense_annotation_base.py")
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)
OBSERVE_MOTIF_RUNTIME = M.motif_runtime
RAW = ("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n##contig=<ID=chrM,length=100>\n"
       "##INFO=<ID=AM_PATHOGENICITY,Number=1,Type=Float,Description=\"AM\">\n"
       "##INFO=<ID=TRANSCRIPT,Number=1,Type=String,Description=\"tx\">\n"
       "##INFO=<ID=PVAR,Number=1,Type=String,Description=\"AA\">\n")
CSQ = '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Gene|Feature|Consequence">\n'
COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
ROWS = ["chr1\t10\t.\tA\tC\t.\tPASS\tAM_PATHOGENICITY=0.8;TRANSCRIPT=ENST1;PVAR=A2T",
        "chrM\t20\t.\tG\tT\t.\tPASS\tAM_PATHOGENICITY=0.9;TRANSCRIPT=ENSTM;PVAR=G3V"]


class AnnotationBaseFixtures(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(prefix="alpha_base_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(temp.cleanup)
        self.folder = Path(temp.name)
        self.cache = self.folder / "cache"
        self.cache.mkdir()
        self.source = self.folder / "original"
        self.scripts = self.source / "scripts"
        self.scripts.mkdir(parents=True)
        (self.scripts / "alphmis_tsv2vcf.py").write_text("# original converter fixture\n")
        (self.scripts / "install_utils.sh").write_text(
            "function basic_vep_annotation() {\n    vep --fork 22 --buffer_size 10000 --use_transcript_ref --total_length --domains\n}\n"
            'function AlphaMissense_anno() {\n    python "$convert_py" "$alphamissense_prescore" "$converted" && \\\n'
            '    basic_vep_annotation \\\n    -i "$converted" -a "$assembly" -v "$conservation_file" && \\\n    :\n}\n')
        (self.scripts / "am_pick_missense_constrained_segments.py").write_text("# unchanged scientific motif fixture\n")
        self.bin = self.folder / "vep_api"
        self.bin.mkdir()
        for name in ("modules/Bio/EnsEMBL/VEP/Runner.pm", "Bio/EnsEMBL/Registry.pm"):
            path = self.bin / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("# pinned VEP API\n")
        self.vep = self.bin / "vep"
        self.vep.write_text("#!/bin/sh\nprintf 'ensembl-vep : 116.3\\n'\n")
        self.vep.chmod(0o755)
        loaded = {"Bio/EnsEMBL/VEP/Runner.pm": str(self.bin / "modules/Bio/EnsEMBL/VEP/Runner.pm"),
                  "Bio/EnsEMBL/Registry.pm": str(self.bin / "Bio/EnsEMBL/Registry.pm")}
        perl = self.bin / "perl"
        perl.write_text("#!" + sys.executable + "\nprint(" + repr(json.dumps({"version": "5.032001", "loaded": loaded})) + ")\n")
        perl.chmod(0o755)
        env = mock.patch.dict(os.environ, {"PATH": str(self.bin) + os.pathsep + os.environ["PATH"],
                                           "PYTHONDONTWRITEBYTECODE": "1"})
        env.start()
        self.addCleanup(env.stop)
        runtime = mock.patch.object(M, "motif_runtime", return_value={"fixture": "same six scientific package versions"})
        runtime.start()
        self.addCleanup(runtime.stop)
        plugin = self.cache / "plugins"
        plugin.mkdir()
        (plugin / "Conservation.pm").write_text("# actual selected plugin fixture\n")
        vep_cache = self.cache / "VEP/homo_sapiens_merged/116_GRCh38"
        vep_cache.mkdir(parents=True)
        (vep_cache / "info.txt").write_text("selected offline merged cache fixture\n")
        self.config = self.folder / "config.json"
        self.cfg = dict(cache_root=str(self.cache), assembly="hg38", vep_installed_version="116", threads=22,
                        alphamissense_vcf=str(self.cache / "raw.vcf.gz"),
                        alphamissense_vep_base_vcf=str(self.cache / "complete.base.vcf.gz"),
                        alphamissense_vep_vcf=str(self.cache / "runtime.vcf.gz"),
                        alphamissense_prescore=str(self.cache / "original.tsv.gz"),
                        ref_genome=str(self.cache / "reference.fa"), conservation_file=str(self.cache / "gerp.bw"),
                        vep_plugins_dir=str(plugin), vep_cache_dir=str(self.cache / "VEP"))
        for name in ("original.tsv.gz", "original.tsv.gz.tbi", "reference.fa", "reference.fa.fai", "gerp.bw"):
            (self.cache / name).write_bytes(b"unchanged scientific resource fixture\n")
        self.write_config()
        for name in ("alphamissense_annotation_base.py", "alphamissense_cache.py"):
            shutil.copyfile(ROOT / "scripts" / name, self.scripts / name)
        # Native pair-publication fixtures isolate the package fingerprint,
        # which is tested below and observed separately in the real environment.
        with (self.scripts / "alphamissense_annotation_base.py").open("a") as stream:
            stream.write('\ndef motif_runtime(): return {"fixture":"same six scientific package versions"}\n')

    def write_config(self):
        self.config.write_text(json.dumps(self.cfg))

    def vcf(self, path, annotated=False, mitochondrial=True, af=False):
        text = RAW
        if annotated:
            text += '##VEP="v116.3" assembly="GRCh38.p14"\n' + CSQ
        if af:
            text += '##INFO=<ID=AF_joint,Number=A,Type=Float,Description="AF">\n'
        rows = ROWS if mitochondrial else ROWS[:1]
        text += COLUMNS
        for row in rows:
            if annotated:
                tx = "ENSTM" if row.startswith("chrM") else "ENST1"
                row += ";CSQ=ENSG|" + tx + "|missense_variant"
            if af:
                row += ";AF_joint=0.1"
            text += row + "\n"
        with Path(path).open("wb") as stream:
            subprocess.run(["bgzip", "-c"], input=text.encode(), stdout=stream, check=True)
        subprocess.run(["tabix", "-f", "-p", "vcf", str(path)], check=True)
        return Path(path)

    def stage(self, mitochondrial=True):
        stage = Path(tempfile.mkdtemp(prefix=".alpha-scientific-stage-", dir=self.cache))
        return self.vcf(stage / "raw.vcf.gz"), self.vcf(stage / "annotated.vcf.gz", True, mitochondrial)

    def publish(self):
        expected = M.capture(self.config, self.source)
        converted, annotation = self.stage()
        M.publish(self.config, self.source, converted, annotation, expected)
        return M.check(self.config, self.source)

    def test_cold_complete_base_is_independent_and_runtime_AF_replacement_keeps_warm_reuse(self):
        before = M.identity(self.cfg["alphamissense_prescore"])
        record = self.publish()
        self.assertEqual(record["counts"], {"chr1": 1, "chrM": 1})
        base, runtime = map(Path, (self.cfg["alphamissense_vep_base_vcf"], self.cfg["alphamissense_vep_vcf"]))
        self.assertEqual(base.read_bytes(), runtime.read_bytes())
        self.assertFalse(base.samefile(runtime))
        binding = M.motif_binding(self.config, self.source)
        self.vcf(runtime, annotated=True, mitochondrial=False, af=True)
        self.assertEqual(M.index_counts(runtime), {"chr1": 1})
        self.assertEqual(M.check(self.config, self.source), record)
        self.assertEqual(M.motif_binding(self.config, self.source), binding)
        self.assertEqual(M.identity(self.cfg["alphamissense_prescore"]), before)

    def test_other_checkout_and_unrelated_shell_edit_reuse_same_scientific_base(self):
        self.publish()
        other = self.folder / "nextflow"
        shutil.copytree(self.source, other)
        with (other / "scripts/install_utils.sh").open("a") as stream:
            stream.write("\nfunction unrelated_installer_fix() {\n    :\n}\n")
        self.assertEqual(M.check(self.config, self.source), M.check(self.config, other))
        self.assertEqual(M.motif_binding(self.config, self.source), M.motif_binding(self.config, other))

    def test_real_scientific_dependencies_and_original_fork_profile_invalidate(self):
        self.publish()
        cases = [Path(self.cfg[key]) for key in ("alphamissense_prescore", "ref_genome", "conservation_file")]
        cases += [self.scripts / "alphmis_tsv2vcf.py", self.bin / "modules/Bio/EnsEMBL/VEP/Runner.pm",
                  Path(self.cfg["vep_plugins_dir"]) / "Conservation.pm",
                  Path(self.cfg["vep_cache_dir"]) / "homo_sapiens_merged/116_GRCh38/info.txt"]
        for path in cases:
            with self.subTest(path=path):
                previous = path.read_bytes()
                path.write_bytes(previous + b"scientific change\n")
                with self.assertRaisesRegex(ValueError, "scientific dependency changed"):
                    M.check(self.config, self.source)
                path.write_bytes(previous)
                # Resource metadata changed legitimately; create a fresh valid
                # tiny build rather than pretending its old timestamp survived.
                self.publish()
        for threads in (14, 21):
            self.cfg["threads"] = threads
            self.write_config()
            with self.assertRaisesRegex(ValueError, "scientific dependency changed"):
                M.check(self.config, self.source)
        self.cfg["threads"] = 22
        self.write_config()
        old = (self.scripts / "install_utils.sh").read_text()
        (self.scripts / "install_utils.sh").write_text(old.replace("10000", "5000"))
        with self.assertRaisesRegex(ValueError, "scientific dependency changed"):
            M.check(self.config, self.source)

    def test_actual_caller_recipe_and_compressed_reference_index_are_bound(self):
        self.publish()
        source = self.scripts / "install_utils.sh"
        original = source.read_text()
        source.write_text(original.replace(' -v "$conservation_file"', ""))
        with self.assertRaisesRegex(ValueError, "scientific dependency changed"):
            M.check(self.config, self.source)
        source.write_text(original)
        self.publish()
        Path(self.cfg["ref_genome"] + ".gzi").write_text("compressed-reference index fixture\n")
        with self.assertRaisesRegex(ValueError, "scientific dependency changed"):
            M.check(self.config, self.source)

    def test_package_versions_and_owned_file_changes_invalidate_only_motif_binding(self):
        self.publish()
        before = M.motif_binding(self.config, self.source)
        with mock.patch.object(M, "motif_runtime", return_value={"fixture": "updated numerical packages"}):
            self.assertNotEqual(M.motif_binding(self.config, self.source), before)
            M.check(self.config, self.source)
        code = self.folder / "package.py"
        code.write_text("original numerical package fixture\n")
        generated = self.folder / "__pycache__/cache.pyc"
        generated.parent.mkdir()
        generated.write_bytes(b"unrelated generated bytecode")
        package = SimpleNamespace(__version__="1.2", __file__=str(code))
        distribution = SimpleNamespace(version="1.2", files=[Path("package.py"), Path("__pycache__/cache.pyc")],
                                       locate_file=lambda name: self.folder / name)
        with mock.patch.object(M.importlib, "import_module", return_value=package), \
             mock.patch.object(M.importlib.metadata, "distribution", return_value=distribution):
            observed = OBSERVE_MOTIF_RUNTIME()
            self.assertEqual(set(observed["packages"]), {"numpy", "pandas", "scipy", "sklearn", "pysam", "joblib"})
            generated.write_bytes(b"changed generated bytecode")
            self.assertEqual(OBSERVE_MOTIF_RUNTIME(), observed)
            code.write_text("changed numerical implementation\n")
            self.assertNotEqual(OBSERVE_MOTIF_RUNTIME(), observed)
            package.__version__ = "1.3"
            self.assertEqual(OBSERVE_MOTIF_RUNTIME()["packages"]["numpy"]["version"], "1.3")

    def test_incomplete_annotation_never_replaces_previous_complete_outputs(self):
        self.publish()
        paths = [Path(self.cfg[k]) for k in ("alphamissense_vcf", "alphamissense_vep_base_vcf", "alphamissense_vep_vcf")]
        before = [p.read_bytes() for p in paths]
        expected = M.capture(self.config, self.source)
        converted, annotation = self.stage(mitochondrial=False)
        with self.assertRaisesRegex(ValueError, "lost or added input records"):
            M.publish(self.config, self.source, converted, annotation, expected)
        self.assertEqual([p.read_bytes() for p in paths], before)

    def test_recover_every_ordered_publication_boundary_without_rerunning_science(self):
        for boundary in ("marker", "converted", "converted_index", "base", "base_index", "canonical", "canonical_index", "receipt"):
            with self.subTest(boundary=boundary):
                expected = M.capture(self.config, self.source)
                converted, annotation = self.stage()
                _, _, paths = M.settings(self.config, self.source)
                targets = {**paths, "marker": M.pending_path(paths["base"]), "receipt": M.receipt_path(paths["base"])}
                targets.update({key + "_index": Path(str(value) + ".tbi") for key, value in paths.items()})
                original_replace = os.replace
                fired = False

                def replace(source, target):
                    nonlocal fired
                    result = original_replace(source, target)
                    if Path(target) == targets[boundary] and not fired:
                        fired = True
                        raise OSError("injected completed-writer interruption")
                    return result

                with mock.patch.object(M.os, "replace", replace):
                    with self.assertRaisesRegex(OSError, "injected"):
                        M.publish(self.config, self.source, converted, annotation, expected)
                self.assertTrue(fired)
                M.recover(self.config, self.source)
                self.assertEqual(M.check(self.config, self.source)["counts"], {"chr1": 1, "chrM": 1})
                self.assertFalse(M.pending_path(paths["base"]).exists())

    def test_verified_native_lineage_adopts_full_input_and_refuses_enriched_candidate(self):
        converted, annotation = self.stage()
        before = M.identity(converted), M.identity(annotation)
        lineage = dict(schema=M.SCHEMA, science=M.capture(self.config, self.source)["science"],
                       counts={"chr1": 1, "chrM": 1}, provenance={"closed_native_run": "fixture"})
        for key, path in (("converted", converted), ("annotation", annotation)):
            lineage[key] = dict(path=str(path), sha256=M.digest(path), index_sha256=M.digest(str(path) + ".tbi"))
        witness = self.folder / "native-lineage.json"
        witness.write_text(json.dumps(lineage))
        M.adopt(self.config, self.source, witness)
        self.assertEqual(M.check(self.config, self.source)["counts"], lineage["counts"])
        self.assertEqual((M.identity(converted), M.identity(annotation)), before)
        enriched = self.vcf(annotation.parent / "enriched.vcf.gz", True, False, True)
        lineage["annotation"] = dict(path=str(enriched), sha256=M.digest(enriched), index_sha256=M.digest(str(enriched) + ".tbi"))
        witness.write_text(json.dumps(lineage))
        with self.assertRaisesRegex(ValueError, "chromosome counts differ"):
            M.adopt(self.config, self.source, witness)

    def test_restore_runtime_copy_preserves_complete_base_and_its_receipt(self):
        self.publish()
        base = Path(self.cfg["alphamissense_vep_base_vcf"])
        before = M.identity(base), M.receipt_path(base).read_bytes()
        Path(self.cfg["alphamissense_vep_vcf"]).unlink()
        M.restore(self.config, self.source)
        self.assertEqual((M.identity(base), M.receipt_path(base).read_bytes()), before)
        self.assertEqual(Path(self.cfg["alphamissense_vep_vcf"]).read_bytes(), base.read_bytes())
        self.assertFalse(Path(self.cfg["alphamissense_vep_vcf"]).samefile(base))

    def test_closed_native_import_adds_only_base_and_preserves_successful_AF_runtime(self):
        converted, annotation = self.stage()
        target = Path(self.cfg["alphamissense_vcf"])
        shutil.copyfile(converted, target)
        shutil.copyfile(str(converted) + ".tbi", str(target) + ".tbi")
        runtime = self.vcf(self.cfg["alphamissense_vep_vcf"], True, False, True)
        proof = self.folder / "closed-native-result.json"
        proof.write_text('{"closed":true,"status":0,"fixture":true}\n')
        science = M.capture(self.config, self.source)["science"]
        lineage = dict(schema=M.SCHEMA, science=science, counts={"chr1": 1, "chrM": 1},
                       provenance={"closed_native_run": "fixture"})
        for key, path in (("converted", target), ("annotation", annotation)):
            lineage[key] = dict(path=str(path), sha256=M.digest(path), index_sha256=M.digest(str(path) + ".tbi"))
        lineage["native_runtime"] = dict(vcf=M.identity(runtime), index=M.identity(str(runtime) + ".tbi"),
            annotation_base_sha256=M.digest(annotation), science_sha256=M.json_digest(science),
            closed_native_proof=dict(path=str(proof), sha256=M.digest(proof)))
        witness = self.folder / "closed-native-base-lineage.json"
        witness.write_text(json.dumps(lineage))
        before = [(M.identity(path), path.read_bytes()) for path in (target, runtime, Path(str(runtime) + ".tbi"))]
        M.adopt(self.config, self.source, witness)
        self.assertEqual([(M.identity(path), path.read_bytes()) for path in (target, runtime, Path(str(runtime) + ".tbi"))], before)
        M.runtime_associated(self.config, self.source)
        self.assertEqual(M.check(self.config, self.source)["counts"], {"chr1": 1, "chrM": 1})
        self.assertEqual(M.index_counts(runtime), {"chr1": 1})
        self.assertGreater(Path(self.cfg["alphamissense_vep_base_vcf"]).stat().st_mtime_ns, runtime.stat().st_mtime_ns)
        self.vcf(runtime, True, False, True)
        with self.assertRaisesRegex(ValueError, "Closed native runtime"):
            M.runtime_associated(self.config, self.source)

    def test_patch_change_and_unknown_legacy_base_cannot_establish_reuse(self):
        _, annotation = self.stage()
        shutil.copyfile(annotation, self.cfg["alphamissense_vep_base_vcf"])
        shutil.copyfile(str(annotation) + ".tbi", self.cfg["alphamissense_vep_base_vcf"] + ".tbi")
        with self.assertRaises(FileNotFoundError):
            M.check(self.config, self.source)
        self.publish()
        self.vep.write_text("#!/bin/sh\nprintf 'ensembl-vep : 116.4\\n'\n")
        with self.assertRaisesRegex(ValueError, "scientific dependency changed"):
            M.check(self.config, self.source)

    def test_native_motif_copy_preserves_history_and_normal_installer_recovers_interruption(self):
        self.publish()
        spec = importlib.util.spec_from_file_location("alpha_base_mapped_fixture", ROOT / "scripts/mapped_artifact_publish.py")
        publisher = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(publisher)
        mapping = self.folder / "selected_mapping.tsv"
        mapping.write_text("selected mapping fixture\n")
        anchors = [self.folder / "hcseeker.tsv", self.folder / "rmc.tsv"]
        for path in anchors:
            path.write_text("original anchor fixture\n")
        historical = self.folder / "preserved-first-native-motif.pkl"
        historical.write_bytes(b"exact preserved first full-input motif fixture\n")
        original_record = dict(mapping_identity=publisher.BIOMART.mapping_identity(mapping, anchors),
                               artifact_bytes=historical.stat().st_size, artifact_sha256=M.digest(historical))
        publisher.BIOMART.atomic_json(publisher.receipt_path(historical), original_record)
        lineage = dict(schema=M.MOTIF_SCHEMA, scientific_input=M.motif_binding(self.config, self.source),
                       provenance={"closed_first_native_motif": "fixture"}, artifact_path=str(historical),
                       receipt_sha256=M.digest(publisher.receipt_path(historical)), **original_record)
        witness = self.folder / "native-motif-lineage.json"
        witness.write_text(json.dumps(lineage))
        before = publisher.file_state(historical), publisher.file_state(publisher.receipt_path(historical))
        for boundary in ("marker", "artifact", "receipt"):
            target = self.cache / ("adopted-" + boundary + ".pkl")
            selected = {"marker": publisher.pending_path(target), "artifact": target,
                        "receipt": publisher.receipt_path(target)}
            real_replace = os.replace
            fired = False

            def replace(source, destination):
                nonlocal fired
                result = real_replace(source, destination)
                if Path(destination) == selected[boundary] and not fired:
                    fired = True
                    raise OSError("native adoption interrupted after " + boundary)
                return result

            with mock.patch.object(os, "replace", replace):
                with self.assertRaisesRegex(OSError, "native adoption interrupted"):
                    M.adopt_motif(self.config, self.source, target, witness, mapping, anchors)
            normal = SimpleNamespace(artifact=target, mapping=mapping, input=anchors, cache_root=self.cache,
                output_dir=self.cache, assembly="hg38", config=self.config,
                producer=self.scripts / "am_pick_missense_constrained_segments.py", kind="motifs",
                annotation_base=Path(self.cfg["alphamissense_vep_base_vcf"]), native_lineage=None)
            self.assertTrue(publisher.recover(normal))
            M.motif_check(self.config, self.source, target)
            self.assertEqual(target.read_bytes(), historical.read_bytes())
            self.assertFalse(target.samefile(historical))
            self.assertEqual((publisher.file_state(historical), publisher.file_state(publisher.receipt_path(historical))), before)
            complete = publisher.file_state(target), publisher.file_state(publisher.receipt_path(target))
            M.adopt_motif(self.config, self.source, target, witness, mapping, anchors)
            self.assertEqual((publisher.file_state(target), publisher.file_state(publisher.receipt_path(target))), complete)


if __name__ == "__main__":
    unittest.main()
