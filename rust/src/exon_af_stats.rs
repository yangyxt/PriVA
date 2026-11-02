//! The one-dimensional NumPy 1.26.4 reductions used by the exon-AF producer.
//!
//! Keep the input sample order. In particular, a sequential sum changes the
//! mean's low bits: NumPy uses eight accumulators in blocks of at most 128.
//! This compatibility profile retains the reference SSE2 build's NaN bits,
//! two-lane extrema and 8192-element ufunc buffers. A different NumPy build
//! can select different NaN payloads; this is not a generic NumPy API.
//! Algorithms adapted from NumPy's loops_utils.h.src, selection.cpp and
//! loops_minmax.dispatch.c.src (v1.26.4); see NUMPY_LICENSE.

fn add(left: f64, right: f64) -> f64 {
    if left.is_nan() {
        left
    } else if right.is_nan() {
        right
    } else {
        left + right
    }
}

// The reference compiler inlines one recursive level. Its direct leaf uses
// r1+r0 and r3+r2; the inlined leaf uses r0+r1 and r2+r3. These additions have
// identical finite results, but their first NaN is observable in the cache.
// Make that order explicit so Rust optimization cannot choose a different NaN.
fn pairwise(values: &[f64], direct: bool) -> f64 {
    if values.len() < 8 {
        return values.iter().fold(-0.0, |sum, value| add(sum, *value));
    }
    if values.len() > 128 {
        let half = (values.len() / 2) & !7;
        let left = pairwise(&values[..half], !direct);
        let right = pairwise(&values[half..], !direct);
        return add(left, right);
    }
    let mut lanes: [f64; 8] = values[..8].try_into().unwrap();
    let end = values.len() & !7;
    for chunk in values[8..end].chunks_exact(8) {
        for (sum, value) in lanes.iter_mut().zip(chunk) {
            *sum = add(*sum, *value);
        }
    }
    let left = if direct {
        add(add(lanes[1], lanes[0]), add(lanes[3], lanes[2]))
    } else {
        add(add(lanes[0], lanes[1]), add(lanes[2], lanes[3]))
    };
    let sum = add(left, add(add(lanes[4], lanes[5]), add(lanes[6], lanes[7])));
    values[end..]
        .iter()
        .fold(sum, |sum, value| add(sum, *value))
}

fn mean(values: &[f64]) -> f64 {
    values
        .chunks(8192)
        .fold(0.0, |sum, chunk| add(pairwise(chunk, true), sum))
        / values.len() as f64
}

fn less(a: f64, b: f64) -> bool {
    a < b || (!a.is_nan() && b.is_nan())
}

fn extreme(left: f64, right: f64, maximum: bool) -> f64 {
    if left.is_nan() || (if maximum { left > right } else { left < right }) {
        left
    } else {
        right
    }
}

fn extreme_block(initial: f64, remaining: &[f64], maximum: bool) -> f64 {
    if remaining.is_empty() {
        return initial;
    }
    let width = 2;
    let mut lanes = [initial; 2];
    let count = remaining.len() / width * width;
    for chunk in remaining[..count].chunks_exact(width) {
        for (lane, value) in lanes.iter_mut().zip(chunk) {
            *lane = extreme(*lane, *value, maximum);
        }
    }
    let mut value = if lanes[..width].iter().any(|v| v.is_nan()) {
        f64::NAN
    } else {
        let mut step = width / 2;
        while step > 0 {
            for i in 0..step {
                lanes[i] = extreme(lanes[i], lanes[i + step], maximum);
            }
            step /= 2;
        }
        lanes[0]
    };
    for &next in &remaining[count..] {
        value = extreme(value, next, maximum);
    }
    value
}

fn extrema(values: &[f64], maximum: bool) -> f64 {
    let first = values.len().min(8192);
    let mut value = extreme_block(values[0], &values[1..first], maximum);
    for chunk in values[first..].chunks(8192) {
        value = extreme_block(value, chunk, maximum);
    }
    value
}

fn remember(pivots: &mut Vec<usize>, pivot: usize, requested: usize) {
    if pivot == requested && pivots.len() == 50 {
        *pivots.last_mut().unwrap() = pivot;
    } else if pivot >= requested && pivots.len() < 50 {
        pivots.push(pivot);
    }
}

fn median_five(values: &mut [f64]) -> usize {
    for (left, right) in [(0, 1), (3, 4), (0, 3), (1, 4), (1, 2)] {
        if less(values[right], values[left]) {
            values.swap(left, right);
        }
    }
    if less(values[3], values[2]) {
        if less(values[3], values[1]) { 1 } else { 3 }
    } else {
        2
    }
}

// Match NumPy's introselect tie traversal, including the NaN chosen for the
// median's final NaN check. Merely sorting equal NaNs can change their bits.
// Reference: numpy/core/src/npysort/selection.cpp (v1.26.4).
fn partition(values: &mut [f64], kth: usize, pivots: &mut Vec<usize>) {
    let (mut low, mut high) = (0, values.len() - 1);
    while let Some(&pivot) = pivots.last() {
        if pivot > kth {
            high = pivot - 1;
            break;
        }
        if pivot == kth {
            return;
        }
        low = pivot + 1;
        pivots.pop();
    }
    if kth - low < 3 {
        for start in low..=kth {
            let mut smallest = start;
            for at in start + 1..=high {
                if less(values[at], values[smallest]) {
                    smallest = at;
                }
            }
            values.swap(start, smallest);
        }
        remember(pivots, kth, kth);
        return;
    }
    if kth == values.len() - 1 {
        let mut greatest = low;
        for at in low + 1..values.len() {
            if !less(values[at], values[greatest]) {
                greatest = at;
            }
        }
        values.swap(kth, greatest);
        return;
    }
    let mut depth = 2 * values.len().ilog2() as i32;
    while low + 1 < high {
        let (mut left, mut right) = (low + 1, high);
        if depth > 0 || right - left < 5 {
            let middle = low + (high - low) / 2;
            for (a, b) in [(middle, high), (low, high), (middle, low)] {
                if less(values[b], values[a]) {
                    values.swap(a, b);
                }
            }
            values.swap(middle, low + 1);
        } else {
            let range = &mut values[left..right];
            let count = range.len() / 5;
            for group in 0..count {
                let chosen = median_five(&mut range[group * 5..group * 5 + 5]);
                range.swap(group, group * 5 + chosen);
            }
            if count > 2 {
                partition(&mut range[..count], count / 2, &mut Vec::new());
            }
            values.swap(low, left + count / 2);
            left -= 1;
            right += 1;
        }
        depth -= 1;
        let pivot = values[low];
        loop {
            left += 1;
            while less(values[left], pivot) {
                left += 1;
            }
            right -= 1;
            while less(pivot, values[right]) {
                right -= 1;
            }
            if right < left {
                break;
            }
            values.swap(left, right);
        }
        values.swap(low, right);
        if right != kth {
            remember(pivots, right, kth);
        }
        if right >= kth {
            high = right - 1;
        }
        if right <= kth {
            low = left;
        }
    }
    if high == low + 1 && less(values[high], values[low]) {
        values.swap(low, high);
    }
    remember(pivots, kth, kth);
}

/// Median, mean, maximum and minimum, in the producer's five-tuple order.
/// Empty sample lists are skipped by the caller, as by `calculate_stats`.
pub(crate) fn calculate(values: &[f64]) -> [f64; 4] {
    assert!(!values.is_empty());
    let mut sorted = values.to_vec();
    let middle = values.len() / 2;
    let start = (values.len() - 1) / 2;
    let mut pivots = Vec::new();
    partition(&mut sorted, start, &mut pivots);
    if middle != start {
        partition(&mut sorted, middle, &mut pivots);
    }
    partition(&mut sorted, values.len() - 1, &mut pivots);
    let median = if sorted.last().unwrap().is_nan() {
        *sorted.last().unwrap()
    } else {
        mean(&sorted[start..=middle])
    };
    [
        median,
        mean(values),
        extrema(values, true),
        extrema(values, false),
    ]
}
