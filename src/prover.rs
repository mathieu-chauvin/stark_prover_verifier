//import the field
use crate::field::FieldElement;
use crate::field::P;

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

fn get_power_cycle(r: FieldElement) -> Vec<FieldElement> {
    let mut o = vec![FieldElement::new(1), r];
    while *o.last().unwrap() != FieldElement::new(1) {
        o.push(*o.last().unwrap() * r);
    }
    o.pop();
    o
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
    let xs = get_power_cycle(g2);
    let last_step_position = xs[((steps-1)*extension_factor) as usize];

    // Generate the computational trace

    let computational_trace = generate_computational_trace(a, b, steps);

    let output = computational_trace.last().unwrap();
    println!("Done generating computational trace");
    
    
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
    fn test_get_power_cycle() {
        let r = FieldElement::new(2);
        assert_eq!(*get_power_cycle(r).last().unwrap()*r, FieldElement::new(1));
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








