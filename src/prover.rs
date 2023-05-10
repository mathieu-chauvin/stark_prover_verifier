//import the field
use crate::field::FieldElement;
use crate::field::P;
use crate::poly::Poly;
/*fn multiset_eq_rec(a: &[u64], b: &[u64], out_prev: &[u64] out: &mut [u64]) {
    out = out_prev * a /b; 
}

fn multiset_eq(a: &[u64], b: &[u64]) -> bool {
    let mut out = [0u64];
    multiset_eq_rec(a, b, [1u64; 2], &mut out);
    out[0] == 1 && out[1] == 1
}*/

fn is_a_power_of_2(x: u64) -> bool {
    x != 0 && (x & (x - 1)) == 0
}

pub const extension_factor: u64 = 8;

fn generate_computational_trace(a: &[FieldElement], b: &[FieldElement], steps:u64) -> Vec<FieldElement> {
    let mut computational_trace = vec![FieldElement::new(1)];
    for i in 0..steps as usize {
        computational_trace.push(
            computational_trace[i] * a[i] / b[i]
        );
    }
    computational_trace

}


pub fn make_proof_multiset(a: &[FieldElement], b: &[FieldElement], steps:u64){
    assert!(steps <= 2u64.pow(32) / extension_factor);
    assert!(is_a_power_of_2(steps));

    let precision = steps * extension_factor;

    // Root of unity such that x^precision=1
    let g2 = FieldElement::new(7).pow((P-1) / precision);

    // Root of unity such that x^steps=1
    let skips = precision / steps;
    let g1= g2.pow(skips);

    // Powers of the higher-order root of unity
    let xs = FieldElement::get_power_cycle(g2);
    let last_step_position = xs[((steps-1)*extension_factor) as usize];

    // Generate the computational trace

    let computational_trace = generate_computational_trace(a, b, steps);

    let output = computational_trace.last().unwrap();
    println!("Done generating computational trace");

    // Interpolate the computational trace into a polynomial P, with each step
    // along a successive power of G1
    let computational_trace_polynomial = Poly::inv_fft(&computational_trace, &g1);
    let p_evaluations = Poly::fft(&computational_trace_polynomial, &g2);
    println!("Converted computational steps into a polynomial and low-degree extended it");

    // construct lagrange interpolation polynomials for A and B sets
    let a_polynomial = Poly::inv_fft(&a.to_vec(), &g1);
    let a_evaluations = Poly::fft(&a_polynomial, &g2);

    let b_polynomial = Poly::inv_fft(&b.to_vec(), &g1);
    let b_evaluations = Poly::fft(&b_polynomial, &g2);

    // access p_evaluations as vec at position [(i+extension_factor)%precision]


    // Create the composed polynomial such that
    // C(P(x), P(g1*x), K(x)) = P(g1*x) - P(x) / A(x) * B(x) mod P
    let c_of_p_evaluations: Vec<_> = (0..precision)
    .map(|i| (p_evaluations[((i+extension_factor)%precision) as usize] -
            p_evaluations[i as usize]
            / a_evaluations[i as usize]
            * b_evaluations[i as usize])
            )
    .collect();
    println!("Computed C(P, K) polynomial");

    // Compute interpolant of ((1, input), (x_atlast_step, output))
    let interpolant =Poly::lagrange_interpolation(&[FieldElement::new(1), last_step_position], &[FieldElement::new(1), FieldElement::new(1)]);
    let i_evaluations = xs.iter().map(|x| Poly::eval(&interpolant, *x)).collect::<Vec<_>>();

    let mut a1 = Poly::new([(FieldElement::new(0)); 256]);
    a1.coeffs[0] = FieldElement::new(P-1);
    a1.coeffs[1] = FieldElement::new(1);
    let mut a2 = Poly::new([(FieldElement::new(0)); 256]);
    a2.coeffs[0] = -last_step_position;
    a2.coeffs[1] = FieldElement::new(1);
    let zeropoly2 = a1.mul_poly(a2);
    let z2_evaluations = xs.iter().map(|x| Poly::eval(&zeropoly2, *x)).collect::<Vec<_>>();
    let inv_z2_evaluations =  FieldElement::multi_inv(&z2_evaluations);

    // B = (P - I) / Z2
    let b_evaluations = p_evaluations.iter()
        .zip(i_evaluations.iter())
        .zip(inv_z2_evaluations.iter())
        .map(|((p, i), invq)| ((*p - *i) * *invq))
        .collect::<Vec<_>>();

    println!("Computed B polynomial");


    
    
}

//tests
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_a_power_of_two(){
        assert_eq!(is_a_power_of_2(8), true);
        assert_eq!(is_a_power_of_2(7), false);
    }

    
    #[test]
    fn test_multiset_computational_trace() {
        let a = [FieldElement::new(1), FieldElement::new(2), FieldElement::new(4), FieldElement::new(3)];
        let b = [FieldElement::new(4), FieldElement::new(2), FieldElement::new(3), FieldElement::new(1)];
        let steps = 4;

        let trace = generate_computational_trace(&a, &b, steps);
        
        assert_eq!(*trace.last().unwrap(), FieldElement::new(1));
    }


}








