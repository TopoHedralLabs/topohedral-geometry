//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
//}}}
//{{{ std imports
//}}}
//{{{ dep imports
//}}}
//--------------------------------------------------------------------------------------------------

mod schema;
mod utils;

pub use schema::{
    Body, BodyId, CurveId, EdgeId, EdgeRef, FaceId, FaceRef, HalfedgeId, HalfedgeRef, Modeller,
    PointId, RegionId, RegionRef, RingId, RingRef, ShellId, ShellRef, SurfaceId, UvCurveId,
    VertexId, VertexRef,
};
