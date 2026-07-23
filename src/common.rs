//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

use core::num;
//{{{ crate imports
//}}}
//{{{ std imports
use std::{ops::Index, ops::IndexMut, vec};
//}}}
//{{{ dep imports
use topohedral_linalg::{DMatrix, DVector, VecType::Col, VectorOps};
//}}}
//--------------------------------------------------------------------------------------------------

pub type Vector = DVector<f64>;
pub type Matrix = DMatrix<f64>;

pub fn zero_vec(n: usize) -> Vector {
    Vector::zeros_vec(n, Col)
}

//{{{ fun: inv_homog
/// Performs the perspective map (inverse of homogeneuos map) from homogeneous coordinates to
/// Euclidean coordinates.
pub fn inv_homog(point_w: &Vector) -> Vector {
    let n = point_w.len();
    let mut point = Vector::from_value_vec(0.0, n - 1, Col);
    let w = point_w[n - 1];
    for i in 0..n - 1 {
        point[i] = point_w[i] / w;
    }
    point
}
//}}}
//{{{ fun: homog
/// Performs the inverse perspective map (homogeneous map) from Euclidean coordinates to
/// Homogeneious coordinates.
pub fn homog(
    point: &Vector,
    weight: f64,
) -> Vector {
    let n = point.len();
    let mut point_w = Vector::from_value_vec(0.0, n + 1, Col);
    for i in 0..n {
        point_w[i] = weight * point[i];
    }
    point_w[n] = weight;
    point_w
}
//}}}
//{{{ struct: IndexHelper
/// This struct provides the index-conversion for going between a 1D index and an N-dimensional index.
///
/// Indexing in this struct is implemented such that the leftmost index varies the fastest.
/// In 2D this is equivalent to column-major ordering. So, given a set of dimensions:
/// $$
/// (n_{1}, n_{2}, ..., n_{N})
/// $$
/// The linear index $j$ of a multi-dimensional index $(i_{1}, i_{2}, ..., i_{N})$ is given by:
/// $$
/// j = i_{1} + n_{1} i_{2} + n_{1}n_{2}i_{3} + ... + n_{1}n_{2}...n_{N-1}i_{N}
/// $$
/// And given a linear index $j$, the multi-dimensional index $(i_{1}, i_{2}, ..., i_{N})$ is:
/// $$
/// todo
/// $$
#[repr(C)]
pub struct IndexHelper<const N: usize> {
    /// Dimensions of the N-dimensional array
    dims: [usize; N],
    /// Helper array for converting a multi-dimensional index into a linear index
    lin_helper: [usize; N],
    /// Helper array for converting a linear index into a multi-dimensional index
    tuple_helper: [usize; N],
}
//}}}
//{{{ impl: IndexHelper
impl<const N: usize> IndexHelper<N> {
    pub fn new(dims: &[usize]) -> Self {
        debug_assert!(dims.len() >= N);

        let mut helper = IndexHelper::<N> {
            dims: [0; N],
            lin_helper: [1; N],
            tuple_helper: [0; N],
        };

        helper.dims.clone_from_slice(&dims[0..N]);
        helper.tuple_helper.clone_from_slice(&dims[0..N]);

        for i in 1..N {
            helper.lin_helper[i] = helper.lin_helper[i - 1] * dims[i - 1];
            helper.tuple_helper[i] = helper.tuple_helper[i] * helper.tuple_helper[i - 1];
        }
        helper
    }

    pub fn lin_index(
        &self,
        indices: &[usize],
    ) -> usize {
        debug_assert!(indices.len() >= N);
        let mut idx = 0usize;
        for i in 0..N {
            idx += indices[i] * self.lin_helper[i];
        }
        idx
    }
    //..............................................................................................

    pub fn tuple_index(
        &self,
        idx: usize,
    ) -> [usize; N] {
        let mut out = [0; N];
        out[0] = idx % self.dims[0];
        for i in 1..N {
            let base = (idx / self.tuple_helper[i]) * self.tuple_helper[i];
            out[i] = (idx - base) / self.tuple_helper[i - 1];
        }
        out
    }
    //.............................................................................................
}
//}}}
//{{{ struct: MultiDimArray
pub struct MultDimArray<T, const MAX_SIZE: usize, const DIM: usize>
where
    T: Clone,
{
    data: [T; MAX_SIZE],
    index_helper: IndexHelper<DIM>,
}
//}}}

impl<T, const MAX_SIZE: usize, const DIM: usize> MultDimArray<T, MAX_SIZE, DIM>
where
    T: Clone,
{
    pub fn new(
        init_value: T,
        dims: &[usize],
    ) -> Self {
        let num_elems: usize = dims.iter().product();
        assert!(dims.len() == DIM);
        assert!(num_elems <= MAX_SIZE);

        Self {
            data: std::array::from_fn(|_| init_value.clone()),
            index_helper: IndexHelper::new(dims),
        }
    }

    pub fn as_slice(&self) -> &[T] {
        &self.data
    }

    pub fn as_slice_mut(&mut self) -> &mut [T] {
        &mut self.data
    }
}

impl<T, const MAX_SIZE: usize, const DIM: usize> Index<&[usize; DIM]>
    for MultDimArray<T, MAX_SIZE, DIM>
where
    T: Clone,
{
    type Output = T;

    fn index(
        &self,
        index: &[usize; DIM],
    ) -> &Self::Output {
        let lin_index = self.index_helper.lin_index(index);
        &self.data[lin_index]
    }
}

impl<T, const MAX_SIZE: usize, const DIM: usize> IndexMut<&[usize; DIM]>
    for MultDimArray<T, MAX_SIZE, DIM>
where
    T: Clone,
{
    fn index_mut(
        &mut self,
        index: &[usize; DIM],
    ) -> &mut Self::Output {
        let lin_index = self.index_helper.lin_index(index);
        &mut self.data[lin_index]
    }
}
