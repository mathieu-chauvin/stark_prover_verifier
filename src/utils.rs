// use sha

use sha2::{Digest, Sha256};

pub fn get_pseudorandom_indices(seed: &[u8], modulus: u64, count: usize, exclude_multiples_of: u64) -> Vec<u64> {
    let mut data = seed.to_vec();
    while data.len() < 8 * count {
        let mut hasher = Sha256::new();
        let size_update = if data.len() > 64 { data.len() - 64 } else { data.len()};
        hasher.update(&data[size_update..]);
        data.extend_from_slice(&hasher.finalize());
    }
    if exclude_multiples_of == 0 {
        (0..count)
            .map(|i| u64::from_be_bytes(data[i * 8..(i + 1) * 8].try_into().unwrap()) % modulus)
            .collect()
    } else {
        let real_modulus = modulus * (exclude_multiples_of - 1) / exclude_multiples_of;
        let o: Vec<u64> = (0..count)
            .map(|i| u64::from_be_bytes(data[i * 8..(i + 1) * 8].try_into().unwrap()) % real_modulus)
            .collect();
        o.iter().map(|&x| x + 1 + x / (exclude_multiples_of - 1)).collect()
    }
}

// tests

#[cfg(test)]

mod tests {
    use super::*;

    #[test]
    fn test_get_pseudorandom_indices(){
        let seed = b"seed";
        let modulus = 100;
        let count = 10;
        let exclude_multiples_of = 42;
        let o = get_pseudorandom_indices(seed, modulus, count, exclude_multiples_of);
        println!("{:?}", o);

        assert_eq!(o.len(), count);
        for i in 0..count {
            assert!(o[i] < modulus);
            assert_ne!(o[i] % exclude_multiples_of, 0);
        }
        
    }



}
