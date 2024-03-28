use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::{Curve, Generator, NonZero, Point, Scalar, SecretScalar};

mod laws {
    use crate::{
        as_raw::AsRaw,
        core::{self, *},
        Generator, NonZero,
    };
    use crate::{Point, Scalar};

    pub trait AlwaysNonZero {}
    impl<T> AlwaysNonZero for NonZero<T> {}

    /// If $A$ and $B$ are valid `Point<E>`, then $A + B$ is a valid `Point<E>`
    ///
    /// For `Point<E>` to be valid it needs to meet two conditions:
    /// 1. It has to be on curve
    /// 2. It has to be free of torsion component
    ///
    /// Sum of two points on curve is always a point on curve by definition, so (1) holds.
    ///
    /// Recall that, generally, any point on elliptic curve can be represented as sum of its
    /// components:
    ///
    /// $$P = p_0 \G + p_1 \T_1 + \dots + p_t \T_t$$
    ///
    /// where $\G$ is a group of large prime order, and $\T_{1,\dots,t}$ are torsion small groups.
    /// Then sum of two points can be represented as:
    ///
    /// $$A + B = (a_0 + b_0) \G + (a_1 + b_1) \T_1 + \dots + (a_t + b_t) \T_t$$
    ///
    /// $A$ and $B$ are valid `Point<E>`, so they are torsion free, which means that
    /// $a_{1,\dots,t} = b_{1,\dots,t} = 0$, so their sum is also torsion free:
    ///
    /// $$A + B = (a_0 + b_0) \G$$
    ///
    /// Therefore, (2) holds.
    #[inline]
    pub fn sum_of_points_is_valid_point<E: Curve>(a: &Point<E>, b: &Point<E>) -> Point<E> {
        let sum = Additive::add(a.as_raw(), b.as_raw());
        // Correctness: refer to doc comment of the function
        Point::from_raw_unchecked(sum)
    }

    /// If $A$ and $B$ are valid `Point<E>`, then $A - B$ are valid `Point<E>`
    ///
    /// Please, refer to [`sum_of_points_is_valid_point`], as the proof is pretty much the same.
    #[inline]
    pub fn sub_of_points_is_valid_point<E: Curve>(a: &Point<E>, b: &Point<E>) -> Point<E> {
        let result = Additive::sub(a.as_raw(), b.as_raw());
        // Correctness: refer to doc comment of the function
        Point::from_raw_unchecked(result)
    }

    /// If $A$ is a valid `Point<E>`, then $A + A$ is a valid `Point<E>`
    ///
    /// The proof is the same as for [`sum_of_points_is_valid_point`], just put `B=A`
    #[inline]
    pub fn double_point_is_valid_point<E: Curve>(x: &Point<E>) -> Point<E> {
        let result = Additive::double(x.as_raw());
        // Correctness: refer to doc comment of the function
        Point::from_raw_unchecked(result)
    }

    /// If $A$ is valid `Point<E>`, then $A + G$ is valid `Point<E>`
    #[inline]
    pub fn sum_of_point_and_generator_is_valid_point<E: Curve>(
        a: &Point<E>,
        g: &Generator<E>,
    ) -> Point<E> {
        sum_of_points_is_valid_point(a, &g.to_point())
    }
    /// If $A$ is valid `Point<E>`, then $G + A$ is valid `Point<E>`
    #[inline]
    pub fn sum_of_generator_and_point_is_valid_point<E: Curve>(
        g: &Generator<E>,
        a: &Point<E>,
    ) -> Point<E> {
        sum_of_points_is_valid_point(&g.to_point(), a)
    }

    /// If $A$ is valid `Point<E>`, then $A - G$ is valid `Point<E>`
    pub fn sub_of_point_and_generator_is_valid_point<E: Curve>(
        a: &Point<E>,
        g: &Generator<E>,
    ) -> Point<E> {
        sub_of_points_is_valid_point(a, &g.to_point())
    }
    /// If $A$ is valid `Point<E>`, then $G - A$ is valid `Point<E>`
    pub fn sub_of_generator_and_point_is_valid_point<E: Curve>(
        g: &Generator<E>,
        a: &Point<E>,
    ) -> Point<E> {
        sub_of_points_is_valid_point(&g.to_point(), a)
    }

    /// If $A$ is a valid `Point<E>`, then $-A$ is a valid `Point<E>`
    ///
    /// [`sub_of_points_is_valid_point`] proves that subtraction of two valid `Point<E>` is a
    /// valid `Point<E>`, so $O - A$ is also a valid `Point<E>`
    #[inline]
    pub fn neg_point_is_valid_point<E: Curve>(a: &Point<E>) -> Point<E> {
        let neg = Additive::negate(a.as_raw());
        // Correctness: refer to doc comment of the function
        Point::from_raw_unchecked(neg)
    }

    /// If $n$ is valid `Scalar<E>` and $A$ is valid `Point<E>`, then $n A$ is a valid `Point<E>`
    ///
    /// For `Point<E>` to be valid it needs to meet two conditions:
    /// 1. It has to be on curve
    /// 2. It has to be free of torsion component
    ///
    /// Point on curve multiplied at any integer is always a point on curve by definition, so
    /// (1) holds.
    ///
    /// Recall that, generally, any point on elliptic curve can be represented as sum of its
    /// components:
    ///
    /// $$P = p_0 \G + p_1 \T_1 + \dots + p_t \T_t$$
    ///
    /// where $\G$ is a group of large prime order, and $\T_{1,\dots,t}$ are torsion small groups.
    /// Then multiplication of point at scalar can be represented as:
    ///
    /// $$nA = n a_0 \G + n a_1 \T_1 + \dots + n a_t \T_t$$
    ///
    /// $A$ is valid `Point<E>`, so it is torsion free, which means that $a_{1,\dots,t} = 0$, so
    /// resulting point is also torsion free:
    ///
    /// $$nA = n a_0 \G$$
    ///
    /// Therefore, (2) holds.
    #[inline]
    pub fn mul_of_scalar_at_point_is_valid_point<E: Curve>(
        n: impl AsRef<Scalar<E>>,
        a: &Point<E>,
    ) -> Point<E> {
        let prod = Multiplicative::mul(n.as_ref().as_raw(), a.as_raw());
        // Correctness: refer to doc comment of the function
        Point::from_raw_unchecked(prod)
    }

    /// Same as [`mul_of_scalar_at_point_is_valid_point`] but flipped arguments
    #[inline]
    pub fn mul_of_point_at_scalar_is_valid_point<E: Curve>(
        a: &Point<E>,
        b: impl AsRef<Scalar<E>>,
    ) -> Point<E> {
        mul_of_scalar_at_point_is_valid_point(b, a)
    }

    /// If $n$ is valid `Scalar<E>`, then $n \G$ is valid `Point<E>`
    ///
    /// Proof is the same as in [`mul_of_scalar_at_point_is_valid_point`] with $A = \G$
    #[inline]
    pub fn mul_of_scalar_at_generator_is_valid_point<E: Curve>(
        n: impl AsRef<Scalar<E>>,
        _g: &Generator<E>,
    ) -> Point<E> {
        let prod = Multiplicative::mul(n.as_ref().as_raw(), &core::CurveGenerator);
        // Correctness: refer to doc comment of the function
        Point::from_raw_unchecked(prod)
    }

    /// Same as [`mul_of_scalar_at_generator_is_valid_point`] but flipped arguments
    #[inline]
    pub fn mul_of_generator_at_scalar_is_valid_point<E: Curve>(
        g: &Generator<E>,
        n: impl AsRef<Scalar<E>>,
    ) -> Point<E> {
        mul_of_scalar_at_generator_is_valid_point(n, g)
    }

    /// If $n$ is valid `NonZero<Scalar<E>>` and $A$ is valid `NonZero<Point<E>>`, then $n A$ is a valid `NonZero<Point<E>>`
    ///
    /// As shown in [`mul_of_scalar_at_point_is_valid_point`], $n A$ is a valid `Point<E>`.
    ///
    /// Since $A$ is free of torsion component and non zero, it has order equal to curve `group_order`,
    /// which means (be definition):
    ///
    /// $$\forall n' < \mathit{group\\_order}: n' A \ne O$$
    ///
    /// As $n$ is valid `Scalar<E>`, it's less than curve `group_order`, therefore $n A \ne O$.
    #[inline]
    pub fn mul_of_nonzero_scalar_at_nonzero_point_is_valid_nonzero_point<E: Curve>(
        n: &(impl AsRef<Scalar<E>> + AlwaysNonZero),
        a: &NonZero<Point<E>>,
    ) -> NonZero<Point<E>> {
        let prod = mul_of_scalar_at_point_is_valid_point(n, a);
        // Correctness: refer to doc comment of the function
        NonZero::new_unchecked(prod)
    }

    /// Same as [`mul_of_nonzero_scalar_at_nonzero_point_is_valid_nonzero_point`] but flipped arguments
    #[inline]
    pub fn mul_of_nonzero_point_at_nonzero_scalar_is_valid_nonzero_point<E: Curve>(
        a: &NonZero<Point<E>>,
        n: &(impl AsRef<Scalar<E>> + AlwaysNonZero),
    ) -> NonZero<Point<E>> {
        mul_of_nonzero_scalar_at_nonzero_point_is_valid_nonzero_point(n, a)
    }

    /// If $n$ is valid `NonZero<Scalar<E>>`, then $n \G$ is valid `NonZero<Point<E>>`
    ///
    /// Proof is the same as in [`mul_of_nonzero_scalar_at_nonzero_point_is_valid_nonzero_point`]
    #[inline]
    pub fn mul_of_nonzero_scalar_at_generator_is_valid_nonzero_point<E: Curve>(
        n: &(impl AsRef<Scalar<E>> + AlwaysNonZero),
        g: &Generator<E>,
    ) -> NonZero<Point<E>> {
        let prod = mul_of_scalar_at_generator_is_valid_point(n, g);
        // Correctness: refer to doc comment of the function
        NonZero::new_unchecked(prod)
    }

    /// Same as [`mul_of_nonzero_scalar_at_generator_is_valid_nonzero_point`] but flipped arguments
    #[inline]
    pub fn mul_of_generator_at_nonzero_scalar_is_valid_nonzero_point<E: Curve>(
        g: &Generator<E>,
        n: &(impl AsRef<Scalar<E>> + AlwaysNonZero),
    ) -> NonZero<Point<E>> {
        mul_of_nonzero_scalar_at_generator_is_valid_nonzero_point(n, g)
    }

    /// If $A$ is valid `NonZero<Point<E>>`, then $-A$ is valid `NonZero<Point<E>>`
    ///
    /// As shown in [`neg_point_is_valid_point`], $-A$ is a valid `Point<E>`.
    ///
    /// Since $A$ is not zero, $-A$ is not zero as well.
    #[inline]
    pub fn neg_nonzero_point_is_nonzero_point<E: Curve>(
        a: &NonZero<Point<E>>,
    ) -> NonZero<Point<E>> {
        let neg = neg_point_is_valid_point(a);
        // Correctness: refer to doc comment of the function
        NonZero::new_unchecked(neg)
    }

    /// If $A$ and $B$ are non-zero scalars mod prime integer $q$, then $A \cdot B \ne 0 \pmod{q}$
    ///
    /// Product of two non-zero integers mod $q$ can be zero if, and only if, $A \cdot B$ divides $q$.
    /// It's not possible as $q$ is prime and $A,B < q$.
    pub fn non_zero_scalar_at_non_zero_scalar_is_non_zero_scalar<E: Curve>(
        a: &(impl AsRef<Scalar<E>> + AlwaysNonZero),
        b: &(impl AsRef<Scalar<E>> + AlwaysNonZero),
    ) -> NonZero<Scalar<E>> {
        let prod = super::scalar::mul(a, b);
        // Correctness: refer to doc commnet of the function
        NonZero::new_unchecked(prod)
    }
}

mod scalar {
    use crate::as_raw::{AsRaw, FromRaw};
    use crate::{core::*, SecretScalar};
    use crate::{NonZero, Scalar};

    #[inline]
    pub fn add<E: Curve>(a: impl AsRef<Scalar<E>>, b: impl AsRef<Scalar<E>>) -> Scalar<E> {
        let sum = Additive::add(a.as_ref().as_raw(), b.as_ref().as_raw());
        Scalar::from_raw(sum)
    }

    #[inline]
    pub fn sub<E: Curve>(a: impl AsRef<Scalar<E>>, b: impl AsRef<Scalar<E>>) -> Scalar<E> {
        let result = Additive::sub(a.as_ref().as_raw(), b.as_ref().as_raw());
        Scalar::from_raw(result)
    }

    #[inline]
    pub fn mul<E: Curve>(a: impl AsRef<Scalar<E>>, b: impl AsRef<Scalar<E>>) -> Scalar<E> {
        let prod = Multiplicative::mul(a.as_ref().as_raw(), b.as_ref().as_raw());
        Scalar::from_raw(prod)
    }

    #[inline]
    pub fn neg<E: Curve>(a: &Scalar<E>) -> Scalar<E> {
        let result = Additive::negate(a.as_raw());
        Scalar::from_raw(result)
    }

    #[inline]
    pub fn neg_nonzero<E: Curve>(a: &NonZero<Scalar<E>>) -> NonZero<Scalar<E>> {
        let neg = neg(a);
        // Correctness: since `a` is not zero, `-a` is not zero by definition
        NonZero::new_unchecked(neg)
    }

    #[inline]
    pub fn neg_nonzero_secret<E: Curve>(a: &NonZero<SecretScalar<E>>) -> NonZero<SecretScalar<E>> {
        let mut a: Scalar<E> = *a.as_ref();
        a *= -Scalar::one();
        // Correctness: since `a` is not zero, `-a` is not zero by definition
        NonZero::new_unchecked(SecretScalar::new(&mut a))
    }
}

macro_rules! impl_binary_ops {
    ($($op:ident ($lhs:ty, $op_fn:ident, $rhs:ty = $out:ty) $impl_fn:path),+,) => {$(
        impl<E: Curve> $op<$rhs> for $lhs {
            type Output = $out;
            #[inline]
            fn $op_fn(self, rhs: $rhs) -> Self::Output {
                $impl_fn(&self, &rhs)
            }
        }
        impl<E: Curve> $op<&$rhs> for $lhs {
            type Output = $out;
            #[inline]
            fn $op_fn(self, rhs: &$rhs) -> Self::Output {
                $impl_fn(&self, rhs)
            }
        }
        impl<E: Curve> $op<$rhs> for &$lhs {
            type Output = $out;
            #[inline]
            fn $op_fn(self, rhs: $rhs) -> Self::Output {
                $impl_fn(self, &rhs)
            }
        }
        impl<E: Curve> $op<&$rhs> for &$lhs {
            type Output = $out;
            #[inline]
            fn $op_fn(self, rhs: &$rhs) -> Self::Output {
                $impl_fn(self, rhs)
            }
        }
    )+};
}

macro_rules! impl_nonzero_ops {
    ($($op:ident ($lhs:ty, $op_fn:ident, $rhs:ty = $out:ty) $impl_fn:path),+,) => {
        impl_binary_ops! {$(
            $op (NonZero<$lhs>, $op_fn, NonZero<$rhs> = $out) $impl_fn,
            $op ($lhs, $op_fn, NonZero<$rhs> = $out) $impl_fn,
            $op (NonZero<$lhs>, $op_fn, $rhs = $out) $impl_fn,
        )+}
    };
}

macro_rules! impl_unary_ops {
    ($($op:ident ($op_fn:ident $ty:ty) $impl_fn:path),*,) => {$(
        impl<E: Curve> $op for $ty {
            type Output = $ty;
            #[inline]
            fn $op_fn(self) -> Self::Output {
                $impl_fn(&self)
            }
        }
        impl<E: Curve> $op for &$ty {
            type Output = $ty;
            #[inline]
            fn $op_fn(self) -> Self::Output {
                $impl_fn(self)
            }
        }
    )*};
}

macro_rules! impl_op_assign {
    ($($ty:ty, $trait:ident, $rhs:ty, $fn:ident, $op:tt),+,) => {$(
        impl<E: Curve> $trait<$rhs> for $ty {
            fn $fn(&mut self, rhs: $rhs) {
                *self = *self $op rhs;
            }
        }
        impl<E: Curve> $trait<&$rhs> for $ty {
            fn $fn(&mut self, rhs: &$rhs) {
                *self = *self $op rhs;
            }
        }
    )+};
}

// Point <> Point, Point <> Scalar, Scalar <> Scalar arithmetic ops
impl_binary_ops! {
    Add (Point<E>, add, Point<E> = Point<E>) laws::sum_of_points_is_valid_point,
    Sub (Point<E>, sub, Point<E> = Point<E>) laws::sub_of_points_is_valid_point,

    Add (Point<E>, add, Generator<E> = Point<E>) laws::sum_of_point_and_generator_is_valid_point,
    Add (Generator<E>, add, Point<E> = Point<E>) laws::sum_of_generator_and_point_is_valid_point,
    Sub (Point<E>, sub, Generator<E> = Point<E>) laws::sub_of_point_and_generator_is_valid_point,
    Sub (Generator<E>, sub, Point<E> = Point<E>) laws::sub_of_generator_and_point_is_valid_point,

    Add (Scalar<E>, add, Scalar<E> = Scalar<E>) scalar::add,
    Sub (Scalar<E>, sub, Scalar<E> = Scalar<E>) scalar::sub,
    Mul (Scalar<E>, mul, Scalar<E> = Scalar<E>) scalar::mul,

    Mul (Point<E>, mul, Scalar<E> = Point<E>) laws::mul_of_point_at_scalar_is_valid_point,
    Mul (Scalar<E>, mul, Point<E> = Point<E>) laws::mul_of_scalar_at_point_is_valid_point,
    Mul (Generator<E>, mul, Scalar<E> = Point<E>) laws::mul_of_generator_at_scalar_is_valid_point,
    Mul (Scalar<E>, mul, Generator<E> = Point<E>) laws::mul_of_scalar_at_generator_is_valid_point,
}

// SecretScalar and NonZero<SecretScalar> multiplication, addition, substruction with Scalar,
// NonZero<Scalar>, Point, and NonZero<Point>
impl_binary_ops! {
    Add (SecretScalar<E>, add, Scalar<E> = Scalar<E>) scalar::add,
    Add (Scalar<E>, add, SecretScalar<E> = Scalar<E>) scalar::add,
    Add (SecretScalar<E>, add, NonZero<Scalar<E>> = Scalar<E>) scalar::add,
    Add (NonZero<Scalar<E>>, add, SecretScalar<E> = Scalar<E>) scalar::add,

    Add (NonZero<SecretScalar<E>>, add, Scalar<E> = Scalar<E>) scalar::add,
    Add (Scalar<E>, add, NonZero<SecretScalar<E>> = Scalar<E>) scalar::add,
    Add (NonZero<SecretScalar<E>>, add, NonZero<Scalar<E>> = Scalar<E>) scalar::add,
    Add (NonZero<Scalar<E>>, add, NonZero<SecretScalar<E>> = Scalar<E>) scalar::add,

    Sub (SecretScalar<E>, sub, Scalar<E> = Scalar<E>) scalar::sub,
    Sub (Scalar<E>, sub, SecretScalar<E> = Scalar<E>) scalar::sub,
    Sub (SecretScalar<E>, sub, NonZero<Scalar<E>> = Scalar<E>) scalar::sub,
    Sub (NonZero<Scalar<E>>, sub, SecretScalar<E> = Scalar<E>) scalar::sub,

    Sub (NonZero<SecretScalar<E>>, sub, Scalar<E> = Scalar<E>) scalar::sub,
    Sub (Scalar<E>, sub, NonZero<SecretScalar<E>> = Scalar<E>) scalar::sub,
    Sub (NonZero<SecretScalar<E>>, sub, NonZero<Scalar<E>> = Scalar<E>) scalar::sub,
    Sub (NonZero<Scalar<E>>, sub, NonZero<SecretScalar<E>> = Scalar<E>) scalar::sub,

    Mul (SecretScalar<E>, mul, Scalar<E> = Scalar<E>) scalar::mul,
    Mul (Scalar<E>, mul, SecretScalar<E> = Scalar<E>) scalar::mul,
    Mul (SecretScalar<E>, mul, NonZero<Scalar<E>> = Scalar<E>) scalar::mul,
    Mul (NonZero<Scalar<E>>, mul, SecretScalar<E> = Scalar<E>) scalar::mul,

    Mul (NonZero<SecretScalar<E>>, mul, Scalar<E> = Scalar<E>) scalar::mul,
    Mul (Scalar<E>, mul, NonZero<SecretScalar<E>> = Scalar<E>) scalar::mul,

    Mul (Point<E>, mul, SecretScalar<E> = Point<E>) laws::mul_of_point_at_scalar_is_valid_point,
    Mul (SecretScalar<E>, mul, Point<E> = Point<E>) laws::mul_of_scalar_at_point_is_valid_point,
    Mul (Generator<E>, mul, SecretScalar<E> = Point<E>) laws::mul_of_generator_at_scalar_is_valid_point,
    Mul (SecretScalar<E>, mul, Generator<E> = Point<E>) laws::mul_of_scalar_at_generator_is_valid_point,
    Mul (NonZero<Point<E>>, mul, SecretScalar<E> = Point<E>) laws::mul_of_point_at_scalar_is_valid_point,
    Mul (SecretScalar<E>, mul, NonZero<Point<E>> = Point<E>) laws::mul_of_scalar_at_point_is_valid_point,

    Mul (Point<E>, mul, NonZero<SecretScalar<E>> = Point<E>) laws::mul_of_point_at_scalar_is_valid_point,
    Mul (NonZero<SecretScalar<E>>, mul, Point<E> = Point<E>) laws::mul_of_scalar_at_point_is_valid_point,
}

// NonZero<Point> <> NonZero<Scalar> arithmetic ops
impl_binary_ops! {
    Mul (NonZero<Point<E>>, mul, NonZero<Scalar<E>> = NonZero<Point<E>>) laws::mul_of_nonzero_point_at_nonzero_scalar_is_valid_nonzero_point,
    Mul (NonZero<Scalar<E>>, mul, NonZero<Point<E>> = NonZero<Point<E>>) laws::mul_of_nonzero_scalar_at_nonzero_point_is_valid_nonzero_point,
    Mul (Generator<E>, mul, NonZero<Scalar<E>> = NonZero<Point<E>>) laws::mul_of_generator_at_nonzero_scalar_is_valid_nonzero_point,
    Mul (NonZero<Scalar<E>>, mul, Generator<E> = NonZero<Point<E>>) laws::mul_of_nonzero_scalar_at_generator_is_valid_nonzero_point,

    Mul (NonZero<Point<E>>, mul, NonZero<SecretScalar<E>> = NonZero<Point<E>>) laws::mul_of_nonzero_point_at_nonzero_scalar_is_valid_nonzero_point,
    Mul (NonZero<SecretScalar<E>>, mul, NonZero<Point<E>> = NonZero<Point<E>>) laws::mul_of_nonzero_scalar_at_nonzero_point_is_valid_nonzero_point,
    Mul (Generator<E>, mul, NonZero<SecretScalar<E>> = NonZero<Point<E>>) laws::mul_of_generator_at_nonzero_scalar_is_valid_nonzero_point,
    Mul (NonZero<SecretScalar<E>>, mul, Generator<E> = NonZero<Point<E>>) laws::mul_of_nonzero_scalar_at_generator_is_valid_nonzero_point,
}

// Point <> NonZero<Point>, Scalar <> NonZero<Scalar>,
// NonZero<Point> <> NonZero<Point>, NonZero<Scalar> <> NonZero<Scalar> arithmetic ops
impl_nonzero_ops! {
    Add (Point<E>, add, Point<E> = Point<E>) laws::sum_of_points_is_valid_point,
    Sub (Point<E>, sub, Point<E> = Point<E>) laws::sub_of_points_is_valid_point,

    Add (Scalar<E>, add, Scalar<E> = Scalar<E>) scalar::add,
    Sub (Scalar<E>, sub, Scalar<E> = Scalar<E>) scalar::sub,

    Add (SecretScalar<E>, add, SecretScalar<E> = Scalar<E>) scalar::add,
    Sub (SecretScalar<E>, sub, SecretScalar<E> = Scalar<E>) scalar::sub,
}

// NonZero<Scalar> * NonZero<Scalar>, Scalar * NonZero<Scalar>, NonZero<Scalar> * Scalar
impl_binary_ops! {
    Mul (NonZero<Scalar<E>>, mul, NonZero<Scalar<E>> = NonZero<Scalar<E>>) laws::non_zero_scalar_at_non_zero_scalar_is_non_zero_scalar,
    Mul (Scalar<E>, mul, NonZero<Scalar<E>> = Scalar<E>) scalar::mul,
    Mul (NonZero<Scalar<E>>, mul, Scalar<E> = Scalar<E>) scalar::mul,

    Mul (NonZero<SecretScalar<E>>, mul, NonZero<SecretScalar<E>> = NonZero<Scalar<E>>) laws::non_zero_scalar_at_non_zero_scalar_is_non_zero_scalar,
    Mul (NonZero<Scalar<E>>, mul, NonZero<SecretScalar<E>> = NonZero<Scalar<E>>) laws::non_zero_scalar_at_non_zero_scalar_is_non_zero_scalar,
    Mul (NonZero<SecretScalar<E>>, mul, NonZero<Scalar<E>> = NonZero<Scalar<E>>) laws::non_zero_scalar_at_non_zero_scalar_is_non_zero_scalar,
}

// Point <> NonZero<Scalar>, NonZero<Point> <> Scalar
impl_binary_ops! {
    Mul (Point<E>, mul, NonZero<Scalar<E>> = Point<E>) laws::mul_of_point_at_scalar_is_valid_point,
    Mul (NonZero<Scalar<E>>, mul, Point<E> = Point<E>) laws::mul_of_scalar_at_point_is_valid_point,
    Mul (NonZero<Point<E>>, mul, Scalar<E> = Point<E>) laws::mul_of_point_at_scalar_is_valid_point,
    Mul (Scalar<E>, mul, NonZero<Point<E>> = Point<E>) laws::mul_of_scalar_at_point_is_valid_point,
}

// -Point, -Scalar, -NonZero<Point>, -NonZero<Scalar>
impl_unary_ops! {
    Neg (neg Point<E>) laws::neg_point_is_valid_point,
    Neg (neg Scalar<E>) scalar::neg,
    Neg (neg NonZero<Point<E>>) laws::neg_nonzero_point_is_nonzero_point,
    Neg (neg NonZero<Scalar<E>>) scalar::neg_nonzero,
    Neg (neg NonZero<SecretScalar<E>>) scalar::neg_nonzero_secret,
}

impl_op_assign! {
    Point<E>, AddAssign, Point<E>, add_assign, +,
    Point<E>, AddAssign, NonZero<Point<E>>, add_assign, +,
    Point<E>, AddAssign, Generator<E>, add_assign, +,

    Point<E>, SubAssign, Point<E>, sub_assign, -,
    Point<E>, SubAssign, NonZero<Point<E>>, sub_assign, -,
    Point<E>, SubAssign, Generator<E>, sub_assign, -,

    Point<E>, MulAssign, Scalar<E>, mul_assign, *,
    Point<E>, MulAssign, NonZero<Scalar<E>>, mul_assign, *,
    Point<E>, MulAssign, SecretScalar<E>, mul_assign, *,
    Point<E>, MulAssign, NonZero<SecretScalar<E>>, mul_assign, *,

    Scalar<E>, AddAssign, Scalar<E>, add_assign, +,
    Scalar<E>, AddAssign, NonZero<Scalar<E>>, add_assign, +,
    Scalar<E>, AddAssign, SecretScalar<E>, add_assign, +,
    Scalar<E>, AddAssign, NonZero<SecretScalar<E>>, add_assign, +,

    Scalar<E>, SubAssign, Scalar<E>, sub_assign, -,
    Scalar<E>, SubAssign, NonZero<Scalar<E>>, sub_assign, -,
    Scalar<E>, SubAssign, SecretScalar<E>, sub_assign, -,
    Scalar<E>, SubAssign, NonZero<SecretScalar<E>>, sub_assign, -,

    Scalar<E>, MulAssign, Scalar<E>, mul_assign, *,
    Scalar<E>, MulAssign, NonZero<Scalar<E>>, mul_assign, *,
    Scalar<E>, MulAssign, SecretScalar<E>, mul_assign, *,
    Scalar<E>, MulAssign, NonZero<SecretScalar<E>>, mul_assign, *,

    NonZero<Point<E>>, MulAssign, NonZero<Scalar<E>>, mul_assign, *,
    NonZero<Point<E>>, MulAssign, NonZero<SecretScalar<E>>, mul_assign, *,
    NonZero<Scalar<E>>, MulAssign, NonZero<Scalar<E>>, mul_assign, *,
    NonZero<Scalar<E>>, MulAssign, NonZero<SecretScalar<E>>, mul_assign, *,
}

impl<E: Curve> Point<E> {
    /// Doubles the point, returns `self + self`
    ///
    /// `point.double()` may be more efficient than `point + point` or `2 * point`
    pub fn double(&self) -> Self {
        laws::double_point_is_valid_point(self)
    }
}

#[cfg(test)]
#[allow(dead_code, clippy::redundant_clone)]
fn ensure_ops_implemented<E: Curve>(
    g: Generator<E>,
    point: Point<E>,
    scalar: Scalar<E>,
    non_zero_point: NonZero<Point<E>>,
    non_zero_scalar: NonZero<Scalar<E>>,
    secret_scalar: SecretScalar<E>,
    non_zero_secret_scalar: NonZero<SecretScalar<E>>,
) {
    macro_rules! assert_binary_ops {
        ($($a:ident $op:tt $b:ident => $out:ty),+,) => {$(
            let _: $out = $a.clone() $op $b.clone();
            let _: $out = &$a $op $b.clone();
            let _: $out = $a.clone() $op &$b;
            let _: $out = &$a $op &$b;

            let _: $out = $b.clone() $op $a.clone();
            let _: $out = &$b $op $a.clone();
            let _: $out = $b.clone() $op &$a;
            let _: $out = &$b $op &$a;
        )+};
    }
    macro_rules! assert_unary_ops {
        ($($op:tt $a:ident => $out:ty),+,) => {$(
            let _: $out = $op $a.clone();
            let _: $out = $op &$a;
        )+};
    }

    macro_rules! assert_op_assign {
        ($($a:ident $op:tt $b:ident);+;) => {{$(
            let mut a = $a.clone();
            a $op $b.clone();
            a $op &$b;
        )+}};
    }

    assert_binary_ops!(
        g * scalar => Point<E>,
        point * scalar => Point<E>,
        g * non_zero_scalar => NonZero<Point<E>>,
        non_zero_point * non_zero_scalar => NonZero<Point<E>>,

        g * secret_scalar => Point<E>,
        point * secret_scalar => Point<E>,
        non_zero_point * secret_scalar => Point<E>,

        g * non_zero_secret_scalar => NonZero<Point<E>>,
        point * non_zero_secret_scalar => Point<E>,
        non_zero_point * non_zero_secret_scalar => NonZero<Point<E>>,

        point + point => Point<E>,
        point + non_zero_point => Point<E>,
        non_zero_point + non_zero_point => Point<E>,

        point - point => Point<E>,
        point - non_zero_point => Point<E>,
        non_zero_point - non_zero_point => Point<E>,

        scalar + scalar => Scalar<E>,
        scalar + non_zero_scalar => Scalar<E>,
        non_zero_scalar + non_zero_scalar => Scalar<E>,

        scalar + secret_scalar => Scalar<E>,
        non_zero_scalar + secret_scalar => Scalar<E>,

        scalar + non_zero_secret_scalar => Scalar<E>,
        non_zero_scalar + non_zero_secret_scalar => Scalar<E>,

        scalar - scalar => Scalar<E>,
        scalar - non_zero_scalar => Scalar<E>,
        non_zero_scalar - non_zero_scalar => Scalar<E>,

        scalar - secret_scalar => Scalar<E>,
        non_zero_scalar - secret_scalar => Scalar<E>,

        scalar - non_zero_secret_scalar => Scalar<E>,
        non_zero_scalar - non_zero_secret_scalar => Scalar<E>,

        scalar * scalar => Scalar<E>,
        scalar * non_zero_scalar => Scalar<E>,
        non_zero_scalar * non_zero_scalar => NonZero<Scalar<E>>,

        scalar * secret_scalar => Scalar<E>,
        non_zero_scalar * secret_scalar => Scalar<E>,

        scalar * non_zero_secret_scalar => Scalar<E>,
        non_zero_scalar * non_zero_secret_scalar => NonZero<Scalar<E>>,

        non_zero_secret_scalar + non_zero_secret_scalar => Scalar<E>,
        non_zero_secret_scalar - non_zero_secret_scalar => Scalar<E>,
        non_zero_secret_scalar * non_zero_secret_scalar => NonZero<Scalar<E>>,
    );

    assert_unary_ops!(
        -point => Point<E>,
        -non_zero_point => NonZero<Point<E>>,
        -scalar => Scalar<E>,
        -non_zero_scalar => NonZero<Scalar<E>>,
        -non_zero_secret_scalar => NonZero<SecretScalar<E>>,
    );

    assert_op_assign!(
        point += point;
        point += non_zero_point;
        point += g;

        point -= point;
        point -= non_zero_point;
        point -= g;

        point *= scalar;
        point *= non_zero_scalar;
        point *= secret_scalar;
        point *= non_zero_scalar;

        non_zero_point *= non_zero_scalar;

        scalar += scalar;
        scalar -= scalar;
        scalar *= scalar;

        scalar += non_zero_scalar;
        scalar -= non_zero_scalar;
        scalar *= non_zero_scalar;

        scalar += secret_scalar;
        scalar -= secret_scalar;
        scalar *= secret_scalar;

        scalar += non_zero_secret_scalar;
        scalar -= non_zero_secret_scalar;
        scalar *= non_zero_secret_scalar;

        non_zero_scalar *= non_zero_scalar;
    );
}
