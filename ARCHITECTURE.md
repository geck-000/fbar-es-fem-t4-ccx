# Architecture

How the pieces connect: the deck converter, the patch series, and the three
user elements with their shared operators and result routines.

![Component map](assets/architecture.png)

The Mermaid source below is authoritative for labels and links; regenerate the
image from it in GitDiagram or mermaid.live after an edit. Click handlers work
in GitDiagram and Mermaid-capable editors; GitHub renders the diagram but
ignores the interactions.

```mermaid
flowchart TD

subgraph group_deck["Deck preparation"]
  node_converter["Deck converter<br/>[fbares.py]"]
  node_outputdeck["F-barES deck"]
end

subgraph group_elements["Element formulation"]
  node_u2["U2 edge element<br/>[e_c3d_u2.f]"]
  node_edgeops["Edge operators<br/>[u2edge.f]"]
  node_u3["U3 volumetric element<br/>[e_c3d_u3.f]"]
  node_u3ops["Volumetric operators<br/>[u3vol.f]"]
  node_u3nl["Finite-strain force<br/>[u3nl.f]"]
  node_u3tangent["Finite-strain tangent<br/>[u3nltan.f]"]
  node_u4["U4 base tetrahedron<br/>[e_c3d_u4.f]"]
  node_results2["Element results<br/>[resultsmech_u2.f]"]
  node_results3["Element results<br/>[resultsmech_u3.f]"]
  node_lock["Node-to-element map lock<br/>[fbar_lock.c]"]
end

subgraph group_solver["CalculiX integration"]
  node_patches["Solver patches"]
end

node_engineer(("Analyst"))
node_input["C3D4 deck"]
node_ccx["CalculiX solver"]
node_solution["Displacements and reactions"]

node_engineer -->|"provides"| node_input
node_input -->|"reads deck"| node_converter
node_converter -->|"writes deck"| node_outputdeck
node_converter -->|"defines connectivity"| node_u2
node_converter -->|"defines support"| node_u3
node_converter -->|"retypes tetrahedra"| node_u4
node_outputdeck -->|"submits deck"| node_ccx
node_patches -->|"enables element path"| node_ccx
node_ccx -->|"dispatches U2"| node_u2
node_ccx -->|"dispatches U3"| node_u3
node_ccx -->|"dispatches U4"| node_u4
node_u2 -->|"uses"| node_edgeops
node_u3 -->|"uses"| node_u3ops
node_u3 -->|"finite-strain force"| node_u3nl
node_u3 -->|"finite-strain tangent"| node_u3tangent
node_u2 -->|"provides results"| node_results2
node_u3 -->|"provides results"| node_results3
node_edgeops -->|"guards map"| node_lock
node_u3ops -->|"guards map"| node_lock
node_ccx -->|"solves for"| node_solution

click node_converter "https://github.com/geck-000/fbar-es-fem-t4-ccx/blob/main/elements_ccx/fbares.py"
click node_u2 "https://github.com/geck-000/fbar-es-fem-t4-ccx/blob/main/elements_ccx/e_c3d_u2.f"
click node_edgeops "https://github.com/geck-000/fbar-es-fem-t4-ccx/blob/main/elements_ccx/u2edge.f"
click node_u3 "https://github.com/geck-000/fbar-es-fem-t4-ccx/blob/main/elements_ccx/e_c3d_u3.f"
click node_u3ops "https://github.com/geck-000/fbar-es-fem-t4-ccx/blob/main/elements_ccx/u3vol.f"
click node_u3nl "https://github.com/geck-000/fbar-es-fem-t4-ccx/blob/main/elements_ccx/u3nl.f"
click node_u3tangent "https://github.com/geck-000/fbar-es-fem-t4-ccx/blob/main/elements_ccx/u3nltan.f"
click node_u4 "https://github.com/geck-000/fbar-es-fem-t4-ccx/blob/main/elements_ccx/e_c3d_u4.f"
click node_results2 "https://github.com/geck-000/fbar-es-fem-t4-ccx/blob/main/elements_ccx/resultsmech_u2.f"
click node_results3 "https://github.com/geck-000/fbar-es-fem-t4-ccx/blob/main/elements_ccx/resultsmech_u3.f"
click node_lock "https://github.com/geck-000/fbar-es-fem-t4-ccx/blob/main/elements_ccx/fbar_lock.c"
click node_patches "https://github.com/geck-000/fbar-es-fem-t4-ccx/tree/main/patches_ccx"

classDef toneNeutral fill:#f8fafc,stroke:#334155,stroke-width:1.5px,color:#0f172a
classDef toneBlue fill:#dbeafe,stroke:#2563eb,stroke-width:1.5px,color:#172554
classDef toneAmber fill:#fef3c7,stroke:#d97706,stroke-width:1.5px,color:#78350f
classDef toneMint fill:#dcfce7,stroke:#16a34a,stroke-width:1.5px,color:#14532d
classDef toneRose fill:#ffe4e6,stroke:#e11d48,stroke-width:1.5px,color:#881337
classDef toneIndigo fill:#e0e7ff,stroke:#4f46e5,stroke-width:1.5px,color:#312e81
classDef toneTeal fill:#ccfbf1,stroke:#0f766e,stroke-width:1.5px,color:#134e4a
class node_converter,node_outputdeck toneBlue
class node_u2,node_edgeops,node_u3,node_u3ops,node_u3nl,node_u3tangent,node_u4,node_results2,node_results3,node_lock toneAmber
class node_patches toneMint
class node_engineer,node_input,node_ccx,node_solution toneIndigo
```
