// FRI commitment scheme implementation

use crate::field::{FieldElement, P};
use crate::merkle_tree::{merkelize,mk_branch};
use crate::poly::Poly;
use crate::utils::{get_pseudorandom_indices};

fn prove_low_degree(
    values: Vec<FieldElement>,
    root_of_unity: FieldElement,
    maxdeg_plus_1: u128,
    avoid_multiples : u32
) -> Vec<Vec<u8>> {
    
    // If the degree is small enough, just return the values
    if maxdeg_plus_1 <= 16 {
        return values.iter().map(|x| x.to_bytes().to_vec()).collect::<Vec<_>>();
    }

    let xs =FieldElement::get_power_cycle(root_of_unity);
    assert_eq!(values.len(), xs.len());

    // Compute the Merkle root of the values
    let m = merkelize(&values.iter().map(|x| x.to_bytes().to_vec()).collect::<Vec<_>>());

    let quarter_len = xs.len() / 4;

    // construct rows
    let mut x_polys = vec![];

    for i in 0..quarter_len {
        let mut xs_poly = vec![];
        let mut ys_poly = vec![];
        for j in 0..4 {
            xs_poly.push(xs[i + quarter_len * j]);
            ys_poly.push(values[i + quarter_len * j]);
        }
        x_polys.push(Poly::lagrange_interpolation(&xs_poly, &ys_poly));
    }

    // get a random x value
    let special_x = FieldElement::new(u64::from_be_bytes(m[1].clone().try_into().unwrap()));

    // construct column by successive evaluations of rows at special_x
    let column = x_polys.iter().map(|p| p.eval(special_x)).collect::<Vec<_>>();

    // Compute the Merkle root of the column
    let m2 = merkelize(&column.iter().map(|x| x.to_bytes().to_vec()).collect::<Vec<_>>());

    let ys = get_pseudorandom_indices(&m2[1], column.len() as u32, 40, avoid_multiples);

    // Compute the positions for the values in the polynomial
    let poly_positions = ys.iter().flat_map(|&y| (1..4).map(|i| y + i*values.len() as u32/4)).collect::<Vec<_>>();

    let ys_branches = ys.iter().map(|&y| mk_branch(&m2, y as usize)).collect::<Vec<Vec<_>>>();

    let positions_branches = poly_positions.iter().map(|&y| mk_branch(&m, y as usize).clone()).collect::<Vec<Vec<_>>>();

    // juxtapose elements of the proof
    let mut proof = vec![];
    for i in 0..4 {
        proof.push(m2[1].clone());
        for j in 0..ys.len() {
            proof.push(ys_branches[j][i].clone());
            proof.push(positions_branches[j][i].clone());
        }
    }

    // Recurse...
    let column_proof = prove_low_degree(
        column,
        root_of_unity.pow(4),
        maxdeg_plus_1 / 4,
        avoid_multiples
    );

    // concatenate the proofs
    proof.extend(column_proof);
    proof
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_proving(){
        // test prove low degree function
        let values = vec![
            FieldElement::new(1),
            FieldElement::new(2),
            FieldElement::new(3),
            FieldElement::new(4),
            FieldElement::new(5),
            FieldElement::new(6),
            FieldElement::new(7),
            FieldElement::new(8),
        ];

        let root_of_unity = FieldElement::new(5);

        let proof = prove_low_degree(values, root_of_unity, 4, 0);

        println!("proof: {:?}", proof);

        assert!(false)
        //assert!(verify_low_degree_proof(&proof, root_of_unity, 4, 0));

    }
}