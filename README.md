# SCDS-PB
**S**patial **C**ache **D**ata **S**tructures with **P**ython **B**ind implementation. The C++ code in this repo is fully templated and I try to incorporate as much modern C++ programming (like perfect-forwarding, rvalue references, conditional compilation and SFINAE), as well as macro programming. This repo is expected to develop high performance CPU/GPU end spatial caching structures for point clouds and graphics primitives, etc. 

Currently supported features:

- CPU end: Static-Multi-Tree (like Quad, Oct-Tree): dynamic insert (deletion not supported yet), fully templated (with visualization). C++/Python API
- CPU end: k-d Tree: tested, dynamic insert (deletion not supported yet). C++/Python API. Even faster than Quad-tree (5-10x).
