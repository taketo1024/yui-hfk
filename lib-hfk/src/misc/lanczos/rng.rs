/// A small deterministic generator for the random start blocks.
pub(super) struct SplitMix64(pub(super) u64);

impl SplitMix64 {
    pub(super) fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let z = (self.0 ^ (self.0 >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        let z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // The reference values of SplitMix64 from seed 0.
    #[test]
    fn split_mix_64() {
        let mut rng = SplitMix64(0);
        assert_eq!(rng.next(), 0xe220_a839_7b1d_cdaf);
        assert_eq!(rng.next(), 0x6e78_9e6a_a1b9_65f4);
    }
}
