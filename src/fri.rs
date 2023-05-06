// FRI commitment scheme implementation

use crate::field::{FieldElement, P};
use crate::merkle_tree::merkelize;
use crate::poly::Poly;

fn prove_low_degree(
    values: Vec<FieldElement>,
    root_of_unity: FieldElement,
    maxdeg_plus_1: u128,
) -> Vec<Vec<FieldElement>> {
    
    if maxdeg_plus_1 <= 16 {
        return vec![values];
    }

    let xs =FieldElement::get_power_cycle(root_of_unity);
    assert_eq!(values.len(), xs.len());

    let m = merkelize(&values.iter().map(|x| x.to_bytes().to_vec()).collect::<Vec<_>>());

    let quarter_len = xs.len() / 4;

    // construct rows
    let x_polys = vec![];

    for i in 0..quarter_len {
        let mut xs_poly = vec![];
        let mut ys_poly = vec![];
        for j in 0..4 {
            xs_poly.push(xs[i + quarter_len * j]);
            ys_poly.push(values[i + quarter_len * j]);
        }
        x_polys.push(Poly::lagrange_interpolation(&xs_poly, &ys_poly));
    }


    let special_x = FieldElement::new(u64::from_be_bytes(m[1].try_into().unwrap()));

    // construct column by successive evaluations of rows at special_x
    let column = x_polys.iter().map(|p| p.eval(special_x)).collect::<Vec<_>>();

    let m2 = merkelize(&column.iter().map(|x| x.to_bytes().to_vec()).collect::<Vec<_>>());

    let ys = get_pseudorandom_indices(&m2[1], column.len(), 40, exclude_multiples_of);

    // Compute the positions for the values in the polynomial
    let poly_positions = ys.iter().flat_map(|&y| (0..4).map(move |j| y + (xs.len() / 4) * j)).collect::<Vec<_>>();

    // This component of the proof, including Merkle branches
    let o = vec![
        m2[1].clone(),
        mk_multi_branch(&m2, &ys),
        mk_multi_branch(&m, &poly_positions),
    ];

    // Recurse...
    let column_proof = prove_low_degree(
        column,
        f.exp(root_of_unity, 4),
        maxdeg_plus_1 / 4,
    );
    vec![o] + column_proof
}
