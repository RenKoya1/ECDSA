use crate::utils::{EllipticCurve, FiniteField, Point};
use num_bigint::{BigUint, RandBigInt};
use rand::{self, Rng};
use sha256::{digest, try_digest};

pub struct ECDSA {
    ec: EllipticCurve,
    g: Point,
    q: BigUint,
}

impl ECDSA {
    // generate Private key and Public Key pair
    pub fn generate_key_pair(&self) -> (BigUint, Point) {
        let priv_key = self.generate_priv_key();
        let pub_key = self.generate_pub_key(&priv_key);
        (priv_key, pub_key)
    }

    pub fn generate_priv_key(&self) -> BigUint {
        self.generate_random_number()
    }

    pub fn generate_pub_key(&self, priv_key: &BigUint) -> Point {
        self.ec.scalar_mul(&self.g, priv_key)
    }

    pub fn generate_random_number(&self) -> BigUint {
      let mut rng = rand::thread_rng();  // インスタンスの作成
      rng.gen_biguint_range(&BigUint::from(1u32), &self.q)
    }

    pub fn sign(&self, hash: &BigUint, priv_key: BigUint, k: BigUint) -> (BigUint, BigUint) {
        assert!(*hash < self.q, "hash must be less than q");
        assert!(priv_key < self.q, "priv_key must be less than q");
        assert!(k < self.q, "rundom number k must be less than q");
        let r_point = self.ec.scalar_mul(&self.g, &k);
        if let Point::Coor(r, _) = r_point {
            let s = FiniteField::mult(&r, &priv_key, &self.q);
            let s = FiniteField::add(&s, &hash, &self.q);
            let k_inv = FiniteField::inv_add(&k, &self.q);
            let s = FiniteField::mult(&s, &k_inv, &self.q);
            return (r, s);
        }
        panic!("r_point should not be Identity")
    }

    pub fn verify(&self, hash: &BigUint, pub_key: Point, signature: &(BigUint, BigUint)) -> bool {
        assert!(*hash < self.q, "hash must be less than q");
        let (r, s) = signature;
        let s_inv = FiniteField::inv_mult(&s, &self.q);
        let u1 = FiniteField::mult(&s_inv, hash, &self.q);
        let u2 = FiniteField::mult(&s_inv, r, &self.q);
        let u1a = self.ec.scalar_mul(&self.g, &u1);
        let u1b = self.ec.scalar_mul(&pub_key, &u2);
        let p = self.ec.add(&u1a, &u1b);
        println!("p: {:?}", p);
        if let Point::Coor(xp, _) = p {
            return xp == *r;
        }

        panic!("p should not be Identity")
    }

    pub fn generate_hash(&self, message: &str) -> BigUint {
        let digest = digest(message);
        let hash_bytes = hex::decode(digest).expect("Decoding failed");
        let hash = BigUint::from_bytes_be(&hash_bytes);
        let hash = hash.modpow(&BigUint::from(1u32), &(&self.q - BigUint::from(1u32)));
        let hash = hash + BigUint::from(1u32);
        hash
    }
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_sign_verify() {
        let ec = EllipticCurve {
            a: BigUint::from(2u32),
            b: BigUint::from(2u32),
            p: BigUint::from(17u32),
        };
        let g = Point::Coor(BigUint::from(5u32), BigUint::from(1u32));
        let q = BigUint::from(19u32);

        let ecdsa = ECDSA { ec, g, q };

        let priv_key = BigUint::from(7u32);
        let pub_key = ecdsa.generate_pub_key(&priv_key);
        //select random number k less than q
        //however r, s must not be 0
        let k = BigUint::from(18u32); 

        let message = "Bob transfer 1 BTC to Alice";
        let hash: BigUint = ecdsa.generate_hash(message);
        let signature = ecdsa.sign(&hash, priv_key, k);
        println!("signature: {:?}", signature);

        let verify_result = ecdsa.verify(&hash, pub_key, &signature);
        assert!(verify_result, "verification failed");
    }
}