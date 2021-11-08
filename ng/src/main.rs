mod mesh;

use std::process::Command;

use self::mesh::MeshMaker;

fn main() -> anyhow::Result<()> {
    // This can be made a bit more contact using `ExitStatus::exit_ok`, once
    // that is stable.
    let status = Command::new("cargo")
        .arg("build")
        .args(["--manifest-path", "model/Cargo.toml"])
        .status()?;
    assert!(status.success());

    // TASK: Read up why those calls are unsafe. Make sure calling them is
    //       sound, and document why that is.
    let model = unsafe {
        let lib = libloading::Library::new("model/target/debug/libmodel.so")?;
        let func: libloading::Symbol<ModelFn> = lib.get(b"model")?;
        func()
    };

    let mut mesh = MeshMaker::new();
    let s = model.cube_size;

    // Define a cube
    let v0 = [-s, -s, -s];
    let v1 = [-s, -s, s];
    let v2 = [-s, s, -s];
    let v3 = [-s, s, s];
    let v4 = [s, -s, -s];
    let v5 = [s, -s, s];
    let v6 = [s, s, -s];
    let v7 = [s, s, s];

    // left
    mesh.triangle([v0, v1, v2]);
    mesh.triangle([v2, v1, v3]);

    // right
    mesh.triangle([v4, v6, v5]);
    mesh.triangle([v6, v7, v5]);

    // front
    mesh.triangle([v0, v4, v1]);
    mesh.triangle([v4, v5, v1]);

    // back
    mesh.triangle([v2, v3, v6]);
    mesh.triangle([v6, v3, v7]);

    // bottom
    mesh.triangle([v0, v2, v6]);
    mesh.triangle([v0, v6, v4]);

    // top
    mesh.triangle([v1, v5, v7]);
    mesh.triangle([v1, v7, v3]);

    let mesh = mesh.make();

    println!("Vertices: {:?}", mesh.vertices().collect::<Vec<_>>());
    println!("Triangles: {:?}", mesh.triangles().collect::<Vec<_>>());

    Ok(())
}

type ModelFn = unsafe extern "C" fn() -> fj::Model;
