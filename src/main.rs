
use stark_prover_verifier::field::FieldElement;
use stark_prover_verifier::fri::{prove_low_degree,verify_low_degree_proof};
use stark_prover_verifier::merkle_tree::merkelize;

pub fn main() {
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

        let root_of_unity = FieldElement::nth_root_of_unity(32);

        println!("Began proving");

        let merkle = merkelize(&values.iter().map(|x| x.to_bytes().to_vec()).collect::<Vec<_>>());

        let proof = prove_low_degree(values, root_of_unity, merkle[1].clone(), 32, 7);

        println!("proof: {:?}", proof);

        assert!(verify_low_degree_proof(&merkle[1], &root_of_unity, &proof, 12, 7));

}
