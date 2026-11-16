use topohedral_geometry::{Body, BodyId, Modeller};

#[test]
fn test_modeller() {
    let mut modeller = Modeller::new();
    let first: BodyId = modeller.add_body();
    let second: BodyId = modeller.add_body();
    assert_ne!(first, second);
    assert_eq!(modeller.body(first).unwrap().id(), Some(first));
    assert_eq!(modeller.body(second).unwrap().id(), Some(second));
    assert_eq!(modeller.body_mut(first).unwrap().id(), Some(first));

    let mut empty = Modeller::new();
    assert!(empty.body(first).is_none());
    assert!(empty.body_mut(first).is_none());
    assert!(Body::new().id().is_none());
}
