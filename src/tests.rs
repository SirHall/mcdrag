use crate::{BoundaryLayerCode, MCDRAG};

#[test]
fn fig_26()
{
    // Refer to figure 26 on pg 57

    let settings = MCDRAG::new(5.7, 5.48, 3.0, 0.5, 1.0, 0.754, 0.0, 1.0, 3.34, BoundaryLayerCode::LT);

    // Mach points
    let machs = [
        0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.925, 0.95, 0.975, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.0, 2.2, 2.5,
        3.0, 3.5, 4.0,
    ];

    // Results from paper
    let paper_results = [
        0.112, 0.111, 0.111, 0.112, 0.113, 0.118, 0.135, 0.158, 0.198, 0.281, 0.313, 0.316, 0.308, 0.297, 0.286, 0.276,
        0.267, 0.258, 0.242, 0.228, 0.209, 0.183, 0.162, 0.145,
    ];

    // Generate the drag coefficient at each point
    let generated = machs.map(|mach| settings.gen(mach));

    // Round to three decimal places
    let generated = generated.map(|cd| (cd * 1_000.0).round() / 1_000.0);

    assert_eq!(generated, paper_results);
}
