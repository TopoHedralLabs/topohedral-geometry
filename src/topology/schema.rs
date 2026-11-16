//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use super::utils::{Index, IndexResource, NULL_IDX, SlotVec, is_null};
use crate::{Curve, Vec3, topology::utils::Handle};
//}}}
//{{{ std imports
//}}}
//{{{ dep imports
//}}}
//--------------------------------------------------------------------------------------------------
#[derive(Clone, Copy, Debug)]
pub struct VertexId(Handle);
#[derive(Clone, Copy, Debug)]
pub struct EdgeId(Handle);
#[derive(Clone, Copy, Debug)]
pub struct HalfedgeId(Handle);
#[derive(Clone, Copy, Debug)]
pub struct RingId(Handle);
#[derive(Clone, Copy, Debug)]
pub struct FaceId(Handle);
#[derive(Clone, Copy, Debug)]
pub struct ShellId(Handle);
#[derive(Clone, Copy, Debug)]
pub struct RegionId(Handle);
#[derive(Clone, Copy, Debug)]
pub struct BodyId(Handle);

pub struct Modeller {
    node_id_resource: IndexResource,
    bodies: SlotVec<Body>,
}

impl Modeller {
    pub fn new() -> Modeller {
        Modeller {
            node_id_resource: IndexResource::new(),
            bodies: SlotVec::new(2),
        }
    }
    pub fn add_body(&mut self) {
        let node_id = self.node_id_resource.peek_index();
        self.bodies.insert(Body::new(node_id));
    }
}

pub struct Body {
    id: Index,
    vertices: SlotVec<Vertex>,
    edges: SlotVec<Edge>,
    hafledges: SlotVec<Halfedge>,
    loops: SlotVec<Ring>,
    faces: SlotVec<Face>,
    shells: SlotVec<Shell>,
    regions: SlotVec<Region>,
    points: SlotVec<Vec3>,
    curves: SlotVec<Box<dyn Curve>>,
    uv_curves: SlotVec<Box<dyn Curve>>,
    surfaces: SlotVec<Box<dyn Curve>>,
}

impl Body {
    pub fn new(node_id: Index) -> Body {
        Body {
            id: node_id,
            vertices: SlotVec::new(100),
            edges: SlotVec::new(100),
            hafledges: SlotVec::new(100),
            loops: SlotVec::new(100),
            faces: SlotVec::new(100),
            shells: SlotVec::new(100),
            regions: SlotVec::new(100),
            points: SlotVec::new(100),
            curves: SlotVec::new(100),
            uv_curves: SlotVec::new(100),
            surfaces: SlotVec::new(100),
        }
    }
}

struct Vertex {
    id: Index,
    point: Index,
    edge_out: Index,
}

impl Default for Vertex {
    fn default() -> Self {
        Self {
            id: NULL_IDX,
            point: NULL_IDX,
            edge_out: NULL_IDX,
        }
    }
}

struct Edge {
    id: Index,
    curve: Index,
    hedge_forward: Index,
}

impl Default for Edge {
    fn default() -> Self {
        Self {
            id: NULL_IDX,
            curve: NULL_IDX,
            hedge_forward: NULL_IDX,
        }
    }
}

struct Halfedge {
    id: Index,
    ring: Index,
    hedge_twin: Index,
    hedge_prev: Index,
    hedge_next: Index,
    vertex_end: Index,
    edge: Index,
    uv_curve: Index,
}

impl Default for Halfedge {
    fn default() -> Self {
        Self {
            id: NULL_IDX,
            ring: NULL_IDX,
            hedge_twin: NULL_IDX,
            hedge_prev: NULL_IDX,
            hedge_next: NULL_IDX,
            vertex_end: NULL_IDX,
            edge: NULL_IDX,
            uv_curve: NULL_IDX,
        }
    }
}

struct Ring {
    id: Index,
    face: Index,
    hedge: Index,
}

impl Default for Ring {
    fn default() -> Self {
        Self {
            id: NULL_IDX,
            face: NULL_IDX,
            hedge: NULL_IDX,
        }
    }
}

struct Face {
    id: Index,
    shell_front: Index,
    shell_back: Index,
    rings: Vec<Index>,
    surface: Index,
    sense: bool,
}

impl Default for Face {
    fn default() -> Self {
        Self {
            id: NULL_IDX,
            shell_front: NULL_IDX,
            shell_back: NULL_IDX,
            rings: Vec::new(),
            surface: NULL_IDX,
            sense: false,
        }
    }
}

struct Shell {
    id: Index,
    region: Index,
    faces: Vec<Index>,
}

impl Default for Shell {
    fn default() -> Self {
        Self {
            id: NULL_IDX,
            region: NULL_IDX,
            faces: Vec::new(),
        }
    }
}

struct Region {
    id: Index,
    shells: Vec<Index>,
    solid: bool,
}

impl Default for Region {
    fn default() -> Self {
        Self {
            id: NULL_IDX,
            shells: Vec::new(),
            solid: false,
        }
    }
}
