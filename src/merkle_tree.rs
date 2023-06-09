// implement a merkle tree

// imporr the hash function
use std::hash::{Hash, Hasher};


use sha2::{Digest, Sha256};

pub fn hash_sha(x: Vec<u8>) -> Vec<u8> {
    Sha256::digest(x).to_vec()
}

// takes a list of leaves and returns a merkle tree of length 2*n
pub fn merkelize(L: &Vec<Vec<u8>>) -> Vec<Vec<u8>> {
    println!("begin merkelize");
    println!("L: {:?}", L);
    let mut nodes: Vec<Vec<u8>> = vec![vec![]; L.len() * 2];
    nodes[L.len()..].clone_from_slice(&L);
    for i in (1..L.len()).rev() {
        //takes nodes at 2*I and 2*I+1 and hashes them together
        println!("i: {}", i);
        nodes[i] = hash_sha(nodes[2 * i].iter().chain(nodes[2 * i + 1].iter()).cloned().collect::<Vec<u8>>());
    }
    nodes
}

pub fn mk_branch(tree: &Vec<Vec<u8>>, mut index: usize) -> Vec<Vec<u8>> {
    index += tree.len() / 2;
    let mut o: Vec<Vec<u8>> = vec![tree[index].clone()];
    while index > 1 {
        o.push(if index & 1 == 1 {
            tree[index - 1].clone()
        } else {
            tree[index + 1].clone()
        });
        index /= 2;
    }
    o
}

pub fn verify_branch(root: &Vec<u8>, index: usize, branch: &Vec<Vec<u8>>) -> bool {
    let mut o = branch.first().unwrap().clone();
    for (i, b) in branch.iter().skip(1).enumerate() {
        o = if index & (1 << i) == 0 {
            hash_sha(o.iter().chain(b.iter()).cloned().collect::<Vec<u8>>())
        } else {
            hash_sha(b.iter().chain(o.iter()).cloned().collect::<Vec<u8>>())
        }
    }
    o == *root
}



// tests

#[cfg(test)]
mod tests {
    use super::*;


    #[test]
    fn test_merkelize() {
        let data = vec![b"a".to_vec(),b"b".to_vec(), b"c".to_vec(), b"d".to_vec()];
        let tree = merkelize(&data);
        print!("{:?}", tree);

        let expected_root = hash_sha([hash_sha([b"a".to_vec(), b"b".to_vec()].concat().to_vec()), hash_sha([b"c".to_vec(), b"d".to_vec()].concat().to_vec())].concat().to_vec());


        //let expected_root = hash_sha(&hash_sha([&b"a".to_vec(),&b"b".to_vec()].concat()).extend_from_slice(&hash_sha(&b"c".to_vec(), &b"d".to_vec())));
        assert_eq!(tree[1], expected_root);
    }

    #[test]
    fn test_mk_branch() {
        let data = vec![b"a".to_vec(),b"b".to_vec(), b"c".to_vec(), b"d".to_vec()];
        let tree = merkelize(&data);
        let branch = mk_branch(&tree, 2);
        let expected_branch = vec![tree[6].clone(),tree[7].clone(), tree[2].clone()];
        assert_eq!(branch, expected_branch);
    }

    #[test]
    fn test_verify_branch() {
        let data = vec![b"a".to_vec(),b"b".to_vec(), b"c".to_vec(), b"d".to_vec()];
        let tree = merkelize(&data);
        let branch = mk_branch(&tree, 2);
        let expected_branch = vec![tree[6].clone(),tree[7].clone(), tree[2].clone()];
        assert_eq!(branch, expected_branch);
        assert_eq!(verify_branch(&tree[1], 2, &branch), true);
    }
}

