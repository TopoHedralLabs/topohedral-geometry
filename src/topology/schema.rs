//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use super::utils::{Handle, SlotVec};
use crate::{Curve, Surface, Vec3};

//}}}
//{{{ std imports
//}}}
//{{{ dep imports
//}}}
//--------------------------------------------------------------------------------------------------

//{{{ struct: VertexId
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct VertexId(Handle);
//}}}
//{{{ struct: EdgeId
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct EdgeId(Handle);
//}}}
//{{{ struct: HalfedgeId
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct HalfedgeId(Handle);
//}}}
//{{{ struct: RingId
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct RingId(Handle);
//}}}
//{{{ struct: FaceId
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct FaceId(Handle);
//}}}
//{{{ struct: ShellId
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct ShellId(Handle);
//}}}
//{{{ struct: RegionId
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct RegionId(Handle);
//}}}
//{{{ struct: BodyId
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct BodyId(Handle);
//}}}
//{{{ struct: PointId
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct PointId(Handle);
//}}}
//{{{ struct: CurveId
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct CurveId(Handle);
//}}}
//{{{ struct: UvCurveId
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct UvCurveId(Handle);
//}}}
//{{{ struct: SurfaceId
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct SurfaceId(Handle);
//}}}

//{{{ struct: Modeller
pub struct Modeller {
    bodies: SlotVec<Body>,
}
//}}}
//{{{ impl: Modeller
impl Modeller {
    pub fn new() -> Modeller {
        Modeller {
            bodies: SlotVec::new(2),
        }
    }
    pub fn add_body(&mut self) -> BodyId {
        let handle = self.bodies.insert(Body::new());
        let id = BodyId(handle);
        self.bodies
            .get_mut(handle)
            .expect("body was just inserted")
            .id = Some(id);
        id
    }

    pub fn body(
        &self,
        id: BodyId,
    ) -> Option<&Body> {
        self.bodies.get(id.0)
    }

    pub fn body_mut(
        &mut self,
        id: BodyId,
    ) -> Option<&mut Body> {
        self.bodies.get_mut(id.0)
    }
}
//}}}

//{{{ struct: Body
/// Owns topology and geometry storage. Entity IDs are local to this body;
/// lookups check the slot and generation, but do not check body identity.
pub struct Body {
    id: Option<BodyId>,
    vertices: SlotVec<VertexData>,
    edges: SlotVec<EdgeData>,
    hafledges: SlotVec<HalfedgeData>,
    loops: SlotVec<RingData>,
    faces: SlotVec<FaceData>,
    shells: SlotVec<ShellData>,
    regions: SlotVec<RegionData>,
    points: SlotVec<Vec3>,
    curves: SlotVec<Box<dyn Curve>>,
    uv_curves: SlotVec<Box<dyn Curve>>,
    surfaces: SlotVec<Box<dyn Surface>>,
}
//}}}
//{{{ impl: Body
impl Body {
    pub fn new() -> Body {
        Body {
            id: None,
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

    pub fn id(&self) -> Option<BodyId> {
        self.id
    }

    pub fn vertex(
        &self,
        id: VertexId,
    ) -> Option<VertexRef<'_>> {
        Some(VertexRef {
            body: self,
            id,
            data: self.vertices.get(id.0)?,
        })
    }

    pub fn edge(
        &self,
        id: EdgeId,
    ) -> Option<EdgeRef<'_>> {
        Some(EdgeRef {
            body: self,
            id,
            data: self.edges.get(id.0)?,
        })
    }

    pub fn halfedge(
        &self,
        id: HalfedgeId,
    ) -> Option<HalfedgeRef<'_>> {
        Some(HalfedgeRef {
            body: self,
            id,
            data: self.hafledges.get(id.0)?,
        })
    }

    pub fn ring(
        &self,
        id: RingId,
    ) -> Option<RingRef<'_>> {
        Some(RingRef {
            body: self,
            id,
            data: self.loops.get(id.0)?,
        })
    }

    pub fn face(
        &self,
        id: FaceId,
    ) -> Option<FaceRef<'_>> {
        Some(FaceRef {
            body: self,
            id,
            data: self.faces.get(id.0)?,
        })
    }

    pub fn shell(
        &self,
        id: ShellId,
    ) -> Option<ShellRef<'_>> {
        Some(ShellRef {
            body: self,
            id,
            data: self.shells.get(id.0)?,
        })
    }

    pub fn region(
        &self,
        id: RegionId,
    ) -> Option<RegionRef<'_>> {
        Some(RegionRef {
            body: self,
            id,
            data: self.regions.get(id.0)?,
        })
    }

    pub fn point(
        &self,
        id: PointId,
    ) -> Option<&Vec3> {
        self.points.get(id.0)
    }

    pub fn curve(
        &self,
        id: CurveId,
    ) -> Option<&dyn Curve> {
        self.curves.get(id.0).map(Box::as_ref)
    }

    pub fn uv_curve(
        &self,
        id: UvCurveId,
    ) -> Option<&dyn Curve> {
        self.uv_curves.get(id.0).map(Box::as_ref)
    }

    pub fn surface(
        &self,
        id: SurfaceId,
    ) -> Option<&dyn Surface> {
        self.surfaces.get(id.0).map(Box::as_ref)
    }
}
//}}}

//{{{ struct: VertexData
#[derive(Default)]
struct VertexData {
    id: Option<VertexId>,
    point: Option<PointId>,
    edge_out: Option<EdgeId>,
}
//}}}
//{{{ struct: VertexRef
/// A vertex borrowed from its owning body.
///
/// The body cannot be replaced or mutably borrowed while a view is still used:
///
/// ```compile_fail
/// use topohedral_geometry::{Body, VertexId};
///
/// fn replace_while_borrowed(body: &mut Body, id: VertexId) {
///     let vertex = body.vertex(id).unwrap();
///     *body = Body::new();
///     let _ = vertex.id();
/// }
/// ```
#[derive(Clone, Copy)]
pub struct VertexRef<'a> {
    body: &'a Body,
    id: VertexId,
    data: &'a VertexData,
}
//}}}
//{{{ impl: VertexRef
impl<'a> VertexRef<'a> {
    pub fn id(&self) -> VertexId {
        self.id
    }

    pub fn body(&self) -> &'a Body {
        self.body
    }

    pub fn point(&self) -> Option<&'a Vec3> {
        self.body.point(self.data.point?)
    }

    pub fn edge(&self) -> Option<EdgeRef<'a>> {
        self.body.edge(self.data.edge_out?)
    }
}
//}}}

//{{{ struct: EdgeData
#[derive(Default)]
struct EdgeData {
    id: Option<EdgeId>,
    curve: Option<CurveId>,
    hedge_forward: Option<HalfedgeId>,
}
//}}}
//{{{ struct: EdgeRef
/// An edge borrowed from its owning body.
#[derive(Clone, Copy)]
pub struct EdgeRef<'a> {
    body: &'a Body,
    id: EdgeId,
    data: &'a EdgeData,
}
//}}}
//{{{ impl: EdgeRef
impl<'a> EdgeRef<'a> {
    pub fn id(&self) -> EdgeId {
        self.id
    }

    pub fn body(&self) -> &'a Body {
        self.body
    }

    pub fn curve(&self) -> Option<&'a dyn Curve> {
        self.body.curve(self.data.curve?)
    }

    pub fn halfedge(&self) -> Option<HalfedgeRef<'a>> {
        self.body.halfedge(self.data.hedge_forward?)
    }
}
//}}}

//{{{ struct: HalfedgeData
#[derive(Default)]
struct HalfedgeData {
    id: Option<HalfedgeId>,
    ring: Option<RingId>,
    hedge_twin: Option<HalfedgeId>,
    hedge_prev: Option<HalfedgeId>,
    hedge_next: Option<HalfedgeId>,
    vertex_end: Option<VertexId>,
    edge: Option<EdgeId>,
    uv_curve: Option<UvCurveId>,
}
//}}}
//{{{ struct: HalfedgeRef
/// A halfedge borrowed from its owning body.
#[derive(Clone, Copy)]
pub struct HalfedgeRef<'a> {
    body: &'a Body,
    id: HalfedgeId,
    data: &'a HalfedgeData,
}
//}}}
//{{{ impl: HalfedgeRef
impl<'a> HalfedgeRef<'a> {
    pub fn id(&self) -> HalfedgeId {
        self.id
    }

    pub fn body(&self) -> &'a Body {
        self.body
    }

    pub fn ring(&self) -> Option<RingRef<'a>> {
        self.body.ring(self.data.ring?)
    }

    pub fn twin(&self) -> Option<HalfedgeRef<'a>> {
        self.body.halfedge(self.data.hedge_twin?)
    }

    pub fn prev(&self) -> Option<HalfedgeRef<'a>> {
        self.body.halfedge(self.data.hedge_prev?)
    }

    pub fn next(&self) -> Option<HalfedgeRef<'a>> {
        self.body.halfedge(self.data.hedge_next?)
    }

    pub fn vertex_end(&self) -> Option<VertexRef<'a>> {
        self.body.vertex(self.data.vertex_end?)
    }

    pub fn edge(&self) -> Option<EdgeRef<'a>> {
        self.body.edge(self.data.edge?)
    }

    pub fn uv_curve(&self) -> Option<&'a dyn Curve> {
        self.body.uv_curve(self.data.uv_curve?)
    }
}
//}}}

//{{{ struct: RingData
#[derive(Default)]
struct RingData {
    id: Option<RingId>,
    face: Option<FaceId>,
    hedge: Option<HalfedgeId>,
}
//}}}
//{{{ struct: RingRef
/// A ring borrowed from its owning body.
#[derive(Clone, Copy)]
pub struct RingRef<'a> {
    body: &'a Body,
    id: RingId,
    data: &'a RingData,
}
//}}}
//{{{ impl: RingRef
impl<'a> RingRef<'a> {
    pub fn id(&self) -> RingId {
        self.id
    }

    pub fn body(&self) -> &'a Body {
        self.body
    }

    pub fn face(&self) -> Option<FaceRef<'a>> {
        self.body.face(self.data.face?)
    }

    pub fn halfedge(&self) -> Option<HalfedgeRef<'a>> {
        self.body.halfedge(self.data.hedge?)
    }
}
//}}}

//{{{ struct: FaceData
#[derive(Default)]
struct FaceData {
    id: Option<FaceId>,
    shell_front: Option<ShellId>,
    shell_back: Option<ShellId>,
    rings: Vec<RingId>,
    surface: Option<SurfaceId>,
    sense: bool,
}
//}}}
//{{{ struct: FaceRef
/// A face borrowed from its owning body.
#[derive(Clone, Copy)]
pub struct FaceRef<'a> {
    body: &'a Body,
    id: FaceId,
    data: &'a FaceData,
}
//}}}
//{{{ impl: FaceRef
impl<'a> FaceRef<'a> {
    pub fn id(&self) -> FaceId {
        self.id
    }

    pub fn body(&self) -> &'a Body {
        self.body
    }

    pub fn shell_front(&self) -> Option<ShellRef<'a>> {
        self.body.shell(self.data.shell_front?)
    }

    pub fn shell_back(&self) -> Option<ShellRef<'a>> {
        self.body.shell(self.data.shell_back?)
    }

    pub fn ring_ids(&self) -> &'a [RingId] {
        &self.data.rings
    }

    /// Visits live linked rings in stored order, skipping stale handles.
    pub fn rings(&self) -> impl Iterator<Item = RingRef<'a>> + 'a + use<'a> {
        let body = self.body;
        self.data.rings.iter().filter_map(move |id| body.ring(*id))
    }

    pub fn surface(&self) -> Option<&'a dyn Surface> {
        self.body.surface(self.data.surface?)
    }

    pub fn sense(&self) -> bool {
        self.data.sense
    }
}
//}}}

//{{{ struct: ShellData
#[derive(Default)]
struct ShellData {
    id: Option<ShellId>,
    region: Option<RegionId>,
    faces: Vec<FaceId>,
}
//}}}
//{{{ struct: ShellRef
/// A shell borrowed from its owning body.
#[derive(Clone, Copy)]
pub struct ShellRef<'a> {
    body: &'a Body,
    id: ShellId,
    data: &'a ShellData,
}
//}}}
//{{{ impl: ShellRef
impl<'a> ShellRef<'a> {
    pub fn id(&self) -> ShellId {
        self.id
    }

    pub fn body(&self) -> &'a Body {
        self.body
    }

    pub fn region(&self) -> Option<RegionRef<'a>> {
        self.body.region(self.data.region?)
    }

    pub fn face_ids(&self) -> &'a [FaceId] {
        &self.data.faces
    }

    /// Visits live linked faces in stored order, skipping stale handles.
    pub fn faces(&self) -> impl Iterator<Item = FaceRef<'a>> + 'a + use<'a> {
        let body = self.body;
        self.data.faces.iter().filter_map(move |id| body.face(*id))
    }
}
//}}}

//{{{ struct: RegionData
#[derive(Default)]
struct RegionData {
    id: Option<RegionId>,
    shells: Vec<ShellId>,
    solid: bool,
}
//}}}
//{{{ struct: RegionRef
/// A region borrowed from its owning body.
#[derive(Clone, Copy)]
pub struct RegionRef<'a> {
    body: &'a Body,
    id: RegionId,
    data: &'a RegionData,
}
//}}}
//{{{ impl: RegionRef
impl<'a> RegionRef<'a> {
    pub fn id(&self) -> RegionId {
        self.id
    }

    pub fn body(&self) -> &'a Body {
        self.body
    }

    pub fn shell_ids(&self) -> &'a [ShellId] {
        &self.data.shells
    }

    /// Visits live linked shells in stored order, skipping stale handles.
    pub fn shells(&self) -> impl Iterator<Item = ShellRef<'a>> + 'a + use<'a> {
        let body = self.body;
        self.data
            .shells
            .iter()
            .filter_map(move |id| body.shell(*id))
    }

    pub fn solid(&self) -> bool {
        self.data.solid
    }
}
//}}}

//{{{ mod: tests
#[cfg(test)]
mod tests {

    use super::*;
    use crate::{Line, Plane};

    struct Fixture {
        body: Body,
        vertex: VertexId,
        edge: EdgeId,
        halfedges: [HalfedgeId; 4],
        ring: RingId,
        face: FaceId,
        shells: [ShellId; 2],
        region: RegionId,
    }

    fn vector(values: [f64; 3]) -> Vec3 {
        Vec3::from_col_slice(&values)
    }

    fn assert_vector_eq(
        actual: &Vec3,
        expected: [f64; 3],
    ) {
        for i in 0..3 {
            approx::assert_abs_diff_eq!(actual[i], expected[i], epsilon = 1.0e-12);
        }
    }

    fn linked_body() -> Fixture {
        let mut body = Body::new();
        let point = PointId(body.points.insert(vector([1.0, 2.0, 3.0])));
        let curve = CurveId(body.curves.insert(Box::new(Line::new(
            vector([4.0, 5.0, 6.0]),
            vector([1.0, 0.0, 0.0]),
        ))));
        let uv_curve = UvCurveId(body.uv_curves.insert(Box::new(Line::new(
            vector([7.0, 8.0, 9.0]),
            vector([0.0, 1.0, 0.0]),
        ))));
        let surface = SurfaceId(body.surfaces.insert(Box::new(Plane::new(
            vector([10.0, 11.0, 12.0]),
            vector([0.0, 0.0, 1.0]),
        ))));

        let vertex = VertexId(body.vertices.insert(VertexData::default()));
        let edge = EdgeId(body.edges.insert(EdgeData::default()));
        let halfedges =
            std::array::from_fn(|_| HalfedgeId(body.hafledges.insert(HalfedgeData::default())));
        let ring = RingId(body.loops.insert(RingData::default()));
        let face = FaceId(body.faces.insert(FaceData::default()));
        let shells = std::array::from_fn(|_| ShellId(body.shells.insert(ShellData::default())));
        let region = RegionId(body.regions.insert(RegionData::default()));

        *body.vertices.get_mut(vertex.0).unwrap() = VertexData {
            point: Some(point),
            edge_out: Some(edge),
            ..Default::default()
        };
        *body.edges.get_mut(edge.0).unwrap() = EdgeData {
            curve: Some(curve),
            hedge_forward: Some(halfedges[0]),
            ..Default::default()
        };
        *body.hafledges.get_mut(halfedges[0].0).unwrap() = HalfedgeData {
            ring: Some(ring),
            hedge_twin: Some(halfedges[1]),
            hedge_prev: Some(halfedges[2]),
            hedge_next: Some(halfedges[3]),
            vertex_end: Some(vertex),
            edge: Some(edge),
            uv_curve: Some(uv_curve),
            ..Default::default()
        };
        *body.loops.get_mut(ring.0).unwrap() = RingData {
            face: Some(face),
            hedge: Some(halfedges[0]),
            ..Default::default()
        };
        *body.faces.get_mut(face.0).unwrap() = FaceData {
            shell_front: Some(shells[0]),
            shell_back: Some(shells[1]),
            rings: vec![ring],
            surface: Some(surface),
            sense: true,
            ..Default::default()
        };
        for shell in shells {
            *body.shells.get_mut(shell.0).unwrap() = ShellData {
                region: Some(region),
                faces: vec![face],
                ..Default::default()
            };
        }
        *body.regions.get_mut(region.0).unwrap() = RegionData {
            shells: vec![shells[1], shells[0]],
            solid: true,
            ..Default::default()
        };

        Fixture {
            body,
            vertex,
            edge,
            halfedges,
            ring,
            face,
            shells,
            region,
        }
    }

    #[test]
    fn views_traverse_topology_and_geometry_in_both_directions() {
        let f = linked_body();
        let vertex = f.body.vertex(f.vertex).unwrap();
        // The returned view borrows the body, not this temporary vertex view.
        let edge = f.body.vertex(f.vertex).unwrap().edge().unwrap();
        let halfedge = edge.halfedge().unwrap();
        let ring = halfedge.ring().unwrap();
        let face = ring.face().unwrap();
        let shell = face.shell_front().unwrap();
        let region = shell.region().unwrap();

        assert_eq!(vertex.id(), f.vertex);
        assert_eq!(edge.id(), f.edge);
        assert_eq!(halfedge.id(), f.halfedges[0]);
        assert_eq!(ring.id(), f.ring);
        assert_eq!(face.id(), f.face);
        assert_eq!(shell.id(), f.shells[0]);
        assert_eq!(region.id(), f.region);
        assert_eq!(halfedge.twin().unwrap().id(), f.halfedges[1]);
        assert_eq!(halfedge.prev().unwrap().id(), f.halfedges[2]);
        assert_eq!(halfedge.next().unwrap().id(), f.halfedges[3]);
        assert_eq!(halfedge.vertex_end().unwrap().id(), f.vertex);
        assert_eq!(halfedge.edge().unwrap().id(), f.edge);
        assert_eq!(ring.halfedge().unwrap().id(), f.halfedges[0]);
        assert_eq!(face.shell_back().unwrap().id(), f.shells[1]);
        assert!(face.sense());
        assert!(region.solid());

        for body in [
            vertex.body(),
            edge.body(),
            halfedge.body(),
            ring.body(),
            face.body(),
            shell.body(),
            region.body(),
        ] {
            assert!(std::ptr::eq(body, &f.body));
        }

        assert_vector_eq(vertex.point().unwrap(), [1.0, 2.0, 3.0]);
        assert_vector_eq(&edge.curve().unwrap().eval_point(2.0), [6.0, 5.0, 6.0]);
        assert_vector_eq(
            &halfedge.uv_curve().unwrap().eval_point(2.0),
            [7.0, 10.0, 9.0],
        );
        assert_vector_eq(
            &face.surface().unwrap().eval_point(0.0, 0.0),
            [10.0, 11.0, 12.0],
        );

        // Iterators and ID slices also outlive the temporary views that created them.
        let rings = f.body.face(f.face).unwrap().rings();
        let faces = f.body.shell(f.shells[0]).unwrap().faces();
        let shells = f.body.region(f.region).unwrap().shells();
        let ring_ids = f.body.face(f.face).unwrap().ring_ids();
        assert_eq!(rings.map(|r| r.id()).collect::<Vec<_>>(), vec![f.ring]);
        assert_eq!(faces.map(|f| f.id()).collect::<Vec<_>>(), vec![f.face]);
        assert_eq!(
            shells.map(|s| s.id()).collect::<Vec<_>>(),
            vec![f.shells[1], f.shells[0]]
        );
        assert_eq!(ring_ids, &[f.ring]);
        assert_eq!(shell.face_ids(), &[f.face]);
        assert_eq!(region.shell_ids(), &[f.shells[1], f.shells[0]]);
    }

    #[test]
    fn stale_links_do_not_resolve_to_reused_slots() {
        let mut f = linked_body();
        f.body.edges.delete(f.edge.0);
        let replacement = EdgeId(f.body.edges.insert(EdgeData::default()));
        assert_ne!(replacement, f.edge);
        assert!(f.body.edge(replacement).is_some());
        assert!(f.body.edge(f.edge).is_none());
        assert!(f.body.vertex(f.vertex).unwrap().edge().is_none());
        assert!(f.body.halfedge(f.halfedges[0]).unwrap().edge().is_none());

        f.body.loops.delete(f.ring.0);
        f.body.loops.insert(RingData::default());
        assert!(f.body.face(f.face).unwrap().rings().next().is_none());
        assert_eq!(f.body.face(f.face).unwrap().ring_ids(), &[f.ring]);

        f.body.faces.delete(f.face.0);
        f.body.faces.insert(FaceData::default());
        assert!(f.body.shell(f.shells[0]).unwrap().faces().next().is_none());

        f.body.shells.delete(f.shells[1].0);
        f.body.shells.insert(ShellData::default());
        assert_eq!(
            f.body
                .region(f.region)
                .unwrap()
                .shells()
                .map(|s| s.id())
                .collect::<Vec<_>>(),
            vec![f.shells[0]]
        );

        let data = f.body.vertices.get_mut(f.vertex.0).unwrap();
        data.edge_out = None;
        let point = data.point.take().unwrap();
        f.body.points.delete(point.0);
        f.body.points.insert(vector([20.0, 21.0, 22.0]));
        assert!(f.body.point(point).is_none());
        assert!(f.body.vertex(f.vertex).unwrap().point().is_none());
        assert!(f.body.vertex(f.vertex).unwrap().edge().is_none());
    }

    #[test]
    fn lookups_into_empty_arenas_return_none() {
        let f = linked_body();
        let empty = Body::new();
        assert!(empty.id().is_none());
        assert!(empty.vertex(f.vertex).is_none());
        assert!(empty.edge(f.edge).is_none());
        assert!(empty.halfedge(f.halfedges[0]).is_none());
        assert!(empty.ring(f.ring).is_none());
        assert!(empty.face(f.face).is_none());
        assert!(empty.shell(f.shells[0]).is_none());
        assert!(empty.region(f.region).is_none());

        let vertex = f.body.vertices.get(f.vertex.0).unwrap();
        let edge = f.body.edges.get(f.edge.0).unwrap();
        let halfedge = f.body.hafledges.get(f.halfedges[0].0).unwrap();
        let face = f.body.faces.get(f.face.0).unwrap();
        assert!(empty.point(vertex.point.unwrap()).is_none());
        assert!(empty.curve(edge.curve.unwrap()).is_none());
        assert!(empty.uv_curve(halfedge.uv_curve.unwrap()).is_none());
        assert!(empty.surface(face.surface.unwrap()).is_none());
    }
}
//}}}
