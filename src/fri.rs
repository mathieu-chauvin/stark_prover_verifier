// FRI commitment scheme implementation

use crate::field::{FieldElement, P};
use crate::merkle_tree::{merkelize,mk_branch, verify_branch};
use crate::poly::Poly;
use crate::utils::{get_pseudorandom_indices};
use core::cmp::min;

#[derive(Debug)]
pub struct FRIComponent {
    root: Vec<u8>,
    values: Vec<Vec<u8>>,
    ys_branches : Vec<Vec<Vec<u8>>>,
    positions_branches : Vec<Vec<Vec<u8>>>,
}

pub fn prove_low_degree(
    values: Vec<FieldElement>,
    root_of_unity: FieldElement,
    merkle_root: Vec<u8>,
    maxdeg_plus_1: u128,
    avoid_multiples : u64
) -> Vec<FRIComponent> {

    println!("Starting FRI proof generation");
    // println variables
    println!("maxdeg_plus_1: {}", maxdeg_plus_1);
    println!("values len: {}", values.len());
    
    let mut fri_component = FRIComponent {
        root: vec![],
        values: vec![],
        ys_branches: vec![],
        positions_branches: vec![],
    };

    // If the degree is small enough, just return the values
    if maxdeg_plus_1 <= 16 {
        fri_component.values = values.iter().map(|x| x.to_bytes().to_vec()).collect::<Vec<_>>();
        return vec![fri_component];
        //return values.iter().map(|x| x.to_bytes().to_vec()).collect::<Vec<_>>();
    }

    let xs =FieldElement::get_power_cycle(root_of_unity);
    assert_eq!(values.len(), xs.len());

    // Compute the Merkle root of the values
    let m = merkelize(&values.iter().map(|x| x.to_bytes().to_vec()).collect::<Vec<_>>());

    let quarter_len = xs.len() / 4;

    println!("Done computing Merkle root");
    println!("m[1]: {:?}", m[1]);

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

    println!("Done constructing rows");
    // get a random x value
    let mut rng_merkle: [u8; 8] = [0; 8];
    rng_merkle.copy_from_slice(&merkle_root[0..8]);

    let special_x: FieldElement = FieldElement::new(u64::from_be_bytes(rng_merkle));

    // construct column by successive evaluations of rows at special_x
    let column = x_polys.iter().map(|p| p.eval(special_x)).collect::<Vec<_>>();

    println!("Done constructing column");

    // Compute the Merkle root of the column
    let m2 = merkelize(&column.iter().map(|x| x.to_bytes().to_vec()).collect::<Vec<_>>());

    let ys = get_pseudorandom_indices(&m2[1], column.len() as u64, 40, avoid_multiples);

    println!("Done computing pseudorandom indices");

    // Compute the positions for the values in the polynomial
    // values are at positions y, y + n/4, y + 2n/4, y + 3n/4, y in ys and n = length of values
    let poly_positions = ys.iter().flat_map(|&y| vec![y, y + quarter_len as u64, y + 2 * quarter_len as u64, y + 3 * quarter_len as u64]).collect::<Vec<_>>();

    let ys_branches = ys.iter().map(|&y| mk_branch(&m2, y as usize)).collect::<Vec<Vec<_>>>();

    let positions_branches = poly_positions.iter().map(|&y| mk_branch(&m, y as usize).clone()).collect::<Vec<Vec<_>>>();

    println!("Done computing branches");

    // juxtapose elements of the proof
    //let mut proof = vec![];
     
    fri_component.root = m2[1].clone();
    fri_component.ys_branches = ys_branches.clone();
    fri_component.positions_branches = positions_branches.clone();

    // Recurse...
    let column_proof = prove_low_degree(
        column,
        root_of_unity.pow(4),
        merkle_root.clone(),
        maxdeg_plus_1 / 4,
        avoid_multiples
    );

    let mut proof = vec![fri_component];

    // concatenate the proofs
    proof.extend(column_proof);
    proof
}

pub fn get_branch_value(
    branch: &Vec<Vec<u8>>, 
) -> Vec<u8> {
    let value = branch[0].clone();
    value
}

 
pub fn verify_low_degree_proof(
    merkle_root: &Vec<u8>, 
    root_of_unity: &FieldElement, 
    proof: &Vec<FRIComponent>, 
    maxdeg_plus_1: usize, 
    exclude_multiples_of: u64
) -> bool {

    let mut root1 = merkle_root.clone();

    //let modulus =FieldElement::new(P);

    // test 
    let mut testval = root_of_unity.clone();
    let mut deg_root = 1;
    while testval != FieldElement::new(1) {
        deg_root *= 2;
        testval = testval * testval;
    }

    let mut root_of_unity = root_of_unity.clone();
    let mut maxdeg = maxdeg_plus_1.clone();

    let quartic_roots_of_unity = vec![
        FieldElement::new(1),
        root_of_unity.pow(deg_root / 4),
        root_of_unity.pow(deg_root / 2),
        root_of_unity.pow(deg_root * 3 / 4)
    ];

    for prf_component in proof.iter().take(proof.len() - 1) {
        let mut root2 = prf_component.root.clone();

        // get special x, a pseudorandom value deciding the column we're going to check

        let mut rng_merkle: [u8; 8] = [0; 8];
        rng_merkle.copy_from_slice(&merkle_root[0..8]);
        let special_x = FieldElement::new(u64::from_be_bytes(rng_merkle));

        // get pseudorandom indices, we test on the column we check
        let ys = get_pseudorandom_indices(
            &root2,
            deg_root / 4,
            40,
            exclude_multiples_of
        );

        // get the positions of the values in the polynomial
        // the four positions, for each y, are y, y + n/4, y + 2n/4, y + 3n/4, where n = roudeg
        // they correspond to positions a same "line" we can interpolate our row polynomials from
        let mut poly_positions = Vec::new();
        for y in &ys {
            for j in 0..4 {
                poly_positions.push(y + (deg_root / 4) * j);
            }
        }

        // Verify Merkle branches for columns and poly positions
        for i in 0..poly_positions.len() {
                //column check
            if !verify_branch(&root1, poly_positions[i] as usize, &prf_component.positions_branches[i]) {
                return false;
            }
        }
        let poly_values = prf_component.positions_branches.iter().map(|y| get_branch_value(&y)).collect::<Vec<_>>();

        for i in 0..ys.len() {
            if !verify_branch(&root2, ys[i] as usize, &prf_component.ys_branches[i]) {
                return false;
            }
        }

        let column_values = prf_component.ys_branches.iter().map(|y| get_branch_value(&y)).collect::<Vec<_>>();

        for i in 0..ys.len(){
            //get x coordinates
            // we have 4 x coordinates for each y coordinate
            let y = ys[i];
            let x1 = root_of_unity.pow(y);
            let xcoord = (0..4).map(|j| (quartic_roots_of_unity[j] * x1)).collect::<Vec<_>>();
           
            let row = poly_values[i * 4..(i + 1) * 4]
                .iter()
                .map(|x| FieldElement::new(u64::from_be_bytes(x[..8].try_into().unwrap())))
                .collect::<Vec<_>>();
           
            //TODO : create a field element from bytes function
            let columnval = FieldElement::new(u64::from_be_bytes(column_values[i][..8].try_into().unwrap()));
            
            // do lagrange interpolation
            let poly = Poly::lagrange_interpolation(&xcoord, &row);

            // check if the point from the column belongs to the polynomial
            let point = poly.eval(special_x);
            if point != columnval {
                return false;
            }

        }

        // update root1 and root_of_unity
        root1 = root2.clone();
        root_of_unity = root_of_unity.pow(4).clone();
        maxdeg /= 4;
        deg_root /= 4;
    
    }
    
    // Verify the direct components of the proof
    let values = proof.last().unwrap().values.clone();
    println!("Verifying degree <= {}", maxdeg_plus_1);
    assert!(maxdeg_plus_1 <= 16);
    
    // Check the Merkle root matches up
    let mtree = merkelize(&values);
    assert_eq!(mtree[1], root1);
    
    // Check the degree of the data
    //TODO : exclude_multiples_of, avoid multiples of the generator
    
    let powers = FieldElement::get_power_cycle(root_of_unity);
    let values_full = values.iter().map(|x| FieldElement::new(u64::from_be_bytes(x[..8].try_into().unwrap()))).collect::<Vec<_>>();
    let max_length = min(maxdeg_plus_1, values_full.len());
    let values_short = values_full[..max_length].to_vec();

    let poly_short = Poly::lagrange_interpolation(
        &powers[..max_length],
        &values_short,
    );

    let poly_full = Poly::lagrange_interpolation(
        &powers,
        &values_full,
    );

    for i in 0..values.len() {
        // if both are equal, poly_full is defined by maxdeg_plus_1 points
        assert_eq!(poly_full.eval(powers[i]), poly_short.eval(powers[i]));
    }
    
    
    println!("FRI proof verified");
    true
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
            FieldElement::new(1),
            FieldElement::new(2),
            FieldElement::new(3),
            FieldElement::new(4),
            FieldElement::new(5),
            FieldElement::new(6),
            FieldElement::new(7),
            FieldElement::new(8),
            FieldElement::new(1),
            FieldElement::new(2),
            FieldElement::new(3),
            FieldElement::new(4),
            FieldElement::new(5),
            FieldElement::new(6),
            FieldElement::new(7),
            FieldElement::new(8),
            FieldElement::new(1),
            FieldElement::new(2),
            FieldElement::new(3),
            FieldElement::new(4),
            FieldElement::new(5),
            FieldElement::new(6),
            FieldElement::new(7),
            FieldElement::new(8),
        ];

        let root_of_unity: FieldElement = FieldElement::new(7);

        //map the values to bytes
        //let values: Vec<[u8; 8]> = values.iter().map(|x| x.to_bytes()).collect();

        let merkle = merkelize(&values.iter().map(|x| x.to_bytes().to_vec()).collect::<Vec<_>>());

        println!("Began proving");

        let proof = prove_low_degree(values, root_of_unity, merkle[1].clone(), 32, 7);

        println!("proof: {:?}", proof);

        //assert!(false)
        assert!(verify_low_degree_proof(&merkle[1], &root_of_unity, &proof, 12, 7));

    }
}