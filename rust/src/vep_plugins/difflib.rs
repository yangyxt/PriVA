//! `"".join(difflib.unified_diff(a.splitlines(True), b.splitlines(True),
//! fromfile=..., tofile=...))` with Python's defaults: three context lines,
//! `lineterm="\n"`, no dates, and `SequenceMatcher(None, a, b)` including its
//! automatic junk heuristic (a line of `b` that occurs more than
//! `len(b) // 100 + 1` times is never a match seed once `b` has 200 lines).
use std::collections::HashMap;

/// `str.splitlines(True)`: every Python line boundary ends a line.
pub(super) fn split_lines(text: &str) -> Vec<&str> {
    let mut lines = Vec::new();
    let mut start = 0;
    let mut chars = text.char_indices().peekable();
    while let Some((index, c)) = chars.next() {
        let end = match c {
            '\r' => match chars.peek() {
                Some((next, '\n')) => {
                    let end = next + 1;
                    chars.next();
                    end
                }
                _ => index + 1,
            },
            '\n' | '\u{b}' | '\u{c}' | '\u{1c}' | '\u{1d}' | '\u{1e}' | '\u{85}' | '\u{2028}'
            | '\u{2029}' => index + c.len_utf8(),
            _ => continue,
        };
        lines.push(&text[start..end]);
        start = end;
    }
    if start < text.len() {
        lines.push(&text[start..]);
    }
    lines
}

type Block = (usize, usize, usize);
type Opcode = (char, usize, usize, usize, usize);

struct Matcher<'a> {
    a: &'a [&'a str],
    b: &'a [&'a str],
    b2j: HashMap<&'a str, Vec<usize>>,
}

impl<'a> Matcher<'a> {
    fn new(a: &'a [&'a str], b: &'a [&'a str]) -> Self {
        let mut b2j: HashMap<&str, Vec<usize>> = HashMap::new();
        for (index, line) in b.iter().enumerate() {
            b2j.entry(*line).or_default().push(index);
        }
        if b.len() >= 200 {
            let limit = b.len() / 100 + 1;
            b2j.retain(|_, indexes| indexes.len() <= limit);
        }
        Self { a, b, b2j }
    }

    fn longest_match(&self, alo: usize, ahi: usize, blo: usize, bhi: usize) -> Block {
        let (mut besti, mut bestj, mut bestsize) = (alo, blo, 0);
        let mut lengths: HashMap<usize, usize> = HashMap::new();
        for i in alo..ahi {
            let mut next = HashMap::new();
            for &j in self.b2j.get(self.a[i]).map_or(&[][..], Vec::as_slice) {
                if j < blo {
                    continue;
                }
                if j >= bhi {
                    break;
                }
                let k = j.checked_sub(1).and_then(|p| lengths.get(&p)).unwrap_or(&0) + 1;
                next.insert(j, k);
                if k > bestsize {
                    (besti, bestj, bestsize) = (i + 1 - k, j + 1 - k, k);
                }
            }
            lengths = next;
        }
        // No line is junk (isjunk is None), so equal neighbours always extend
        // the match, popular lines included.
        while besti > alo && bestj > blo && self.a[besti - 1] == self.b[bestj - 1] {
            (besti, bestj, bestsize) = (besti - 1, bestj - 1, bestsize + 1);
        }
        while besti + bestsize < ahi
            && bestj + bestsize < bhi
            && self.a[besti + bestsize] == self.b[bestj + bestsize]
        {
            bestsize += 1;
        }
        (besti, bestj, bestsize)
    }

    fn matching_blocks(&self) -> Vec<Block> {
        let (la, lb) = (self.a.len(), self.b.len());
        let mut queue = vec![(0, la, 0, lb)];
        let mut blocks = Vec::new();
        while let Some((alo, ahi, blo, bhi)) = queue.pop() {
            let (i, j, k) = self.longest_match(alo, ahi, blo, bhi);
            if k > 0 {
                blocks.push((i, j, k));
                if alo < i && blo < j {
                    queue.push((alo, i, blo, j));
                }
                if i + k < ahi && j + k < bhi {
                    queue.push((i + k, ahi, j + k, bhi));
                }
            }
        }
        blocks.sort_unstable();
        let (mut i1, mut j1, mut k1) = (0, 0, 0);
        let mut joined = Vec::new();
        for (i2, j2, k2) in blocks {
            if i1 + k1 == i2 && j1 + k1 == j2 {
                k1 += k2;
            } else {
                if k1 > 0 {
                    joined.push((i1, j1, k1));
                }
                (i1, j1, k1) = (i2, j2, k2);
            }
        }
        if k1 > 0 {
            joined.push((i1, j1, k1));
        }
        joined.push((la, lb, 0));
        joined
    }

    fn opcodes(&self) -> Vec<Opcode> {
        let (mut i, mut j) = (0, 0);
        let mut codes = Vec::new();
        for (ai, bj, size) in self.matching_blocks() {
            let tag = match (i < ai, j < bj) {
                (true, true) => Some('r'),
                (true, false) => Some('d'),
                (false, true) => Some('i'),
                (false, false) => None,
            };
            if let Some(tag) = tag {
                codes.push((tag, i, ai, j, bj));
            }
            (i, j) = (ai + size, bj + size);
            if size > 0 {
                codes.push(('e', ai, i, bj, j));
            }
        }
        codes
    }

    fn grouped_opcodes(&self, n: usize) -> Vec<Vec<Opcode>> {
        let mut codes = self.opcodes();
        if codes.is_empty() {
            codes.push(('e', 0, 1, 0, 1));
        }
        if let Some(first) = codes.first_mut()
            && first.0 == 'e'
        {
            let (tag, i1, i2, j1, j2) = *first;
            *first = (tag, i1.max(i2.saturating_sub(n)), i2, j1.max(j2.saturating_sub(n)), j2);
        }
        if let Some(last) = codes.last_mut()
            && last.0 == 'e'
        {
            let (tag, i1, i2, j1, j2) = *last;
            *last = (tag, i1, i2.min(i1 + n), j1, j2.min(j1 + n));
        }
        let mut groups = Vec::new();
        let mut group = Vec::new();
        for (tag, mut i1, i2, mut j1, j2) in codes {
            if tag == 'e' && i2 - i1 > 2 * n {
                group.push((tag, i1, i2.min(i1 + n), j1, j2.min(j1 + n)));
                groups.push(std::mem::take(&mut group));
                (i1, j1) = (i1.max(i2.saturating_sub(n)), j1.max(j2.saturating_sub(n)));
            }
            group.push((tag, i1, i2, j1, j2));
        }
        if !group.is_empty() && !(group.len() == 1 && group[0].0 == 'e') {
            groups.push(group);
        }
        groups
    }
}

/// `difflib._format_range_unified`
fn range(start: usize, stop: usize) -> String {
    let length = stop - start;
    match length {
        1 => format!("{}", start + 1),
        0 => format!("{start},0"),
        _ => format!("{},{length}", start + 1),
    }
}

pub(super) fn unified_diff(a: &str, b: &str, fromfile: &str, tofile: &str) -> String {
    let (a, b) = (split_lines(a), split_lines(b));
    let mut out = String::new();
    for (index, group) in Matcher::new(&a, &b).grouped_opcodes(3).iter().enumerate() {
        if index == 0 {
            out.push_str(&format!("--- {fromfile}\n+++ {tofile}\n"));
        }
        let (first, last) = (group[0], group[group.len() - 1]);
        out.push_str(&format!(
            "@@ -{} +{} @@\n",
            range(first.1, last.2),
            range(first.3, last.4)
        ));
        for &(tag, i1, i2, j1, j2) in group {
            if tag == 'e' {
                a[i1..i2].iter().for_each(|line| out.push_str(&format!(" {line}")));
                continue;
            }
            if tag == 'r' || tag == 'd' {
                a[i1..i2].iter().for_each(|line| out.push_str(&format!("-{line}")));
            }
            if tag == 'r' || tag == 'i' {
                b[j1..j2].iter().for_each(|line| out.push_str(&format!("+{line}")));
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sha256::sha256_hex;

    // Expected texts are CPython 3.12 difflib output for the same inputs.
    #[test]
    fn hunks_merge_across_six_equal_lines_and_keep_a_missing_final_newline() {
        let a = "one\ntwo\nthree\nfour\nfive\nsix\nseven\neight\nnine\nten\n";
        let b = "one\ntwo\nTHREE\nfour\nfive\nsix\nseven\neight\nnine\nTEN\neleven";
        assert_eq!(
            unified_diff(a, b, "a", "b"),
            "--- a\n+++ b\n@@ -1,10 +1,11 @@\n one\n two\n-three\n+THREE\n four\n five\n six\n \
             seven\n eight\n nine\n-ten\n+TEN\n+eleven"
        );
        assert_eq!(unified_diff("x\ny", "x\ny", "a", "b"), "");
        assert_eq!(unified_diff("", "a\n", "a", "b"), "--- a\n+++ b\n@@ -0,0 +1 @@\n+a\n");
    }

    #[test]
    fn popular_lines_are_not_match_seeds_once_b_has_200_lines() {
        let a: String = (0..240)
            .map(|i| if i % 3 == 0 { "}\n".to_owned() } else { format!("line {i}\n") })
            .collect();
        let b = a.replace("line 100\n", "changed 100\n").replace("line 5\n", "");
        let out = unified_diff(&a, &b, "a", "b");
        assert_eq!(out.len(), 159);
        assert_eq!(
            sha256_hex(out.as_bytes()),
            "b4eb2b477f34fd7d38d22aaa110a589c8888bbe4b81e33e697771f3f361ed244"
        );
    }

    #[test]
    fn the_utr_reporting_patch_matches_the_python_installation() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let read = |path: &str| std::fs::read_to_string(root.join(path)).unwrap();
        let out = unified_diff(
            &read("tests/fixtures/UTRAnnotator116.original.pm"),
            &read("data/vep/UTRAnnotator116_all_effects.pm"),
            "UTRAnnotator.original.pm",
            "UTRAnnotator.pm",
        );
        assert_eq!(out.len(), 4576);
        assert_eq!(
            sha256_hex(out.as_bytes()),
            "f0537be2f1de37d0528014188cbe432558d9409c2018e7708ea7b86d98a5a2c5"
        );
    }

    #[test]
    fn every_python_line_boundary_splits() {
        let a = "a\rb\u{b}c\u{85}d\u{2028}e\r\nf";
        let b = "a\rB\u{b}c\u{85}d\u{2028}e\r\nf";
        assert_eq!(
            unified_diff(a, b, "a", "b"),
            "--- a\n+++ b\n@@ -1,5 +1,5 @@\n a\r-b\u{b}+B\u{b} c\u{85} d\u{2028} e\r\n"
        );
    }
}
