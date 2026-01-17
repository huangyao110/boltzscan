import py3Dmol
from pathlib import Path
from typing import List, Literal, Optional, Tuple

def show_multiple_structures_grid(
    pdb_paths: List[str],
    grid_layout: Optional[Tuple[int, int]] = None,
    show_sidechains: bool = False,
    show_mainchains: bool = False,
    color_scheme: Literal["plddt", "rainbow", "chain"] = "plddt"
) -> py3Dmol.view:
    """
    使用 py3Dmol 以网格布局并排显示多个 PDB 或 CIF 结构文件。

    Args:
        pdb_paths (List[str]): 包含 PDB 或 CIF 文件路径的列表。
        grid_layout (Optional[Tuple[int, int]], optional): 
            一个元组 (rows, cols) 定义网格布局。
            如果为 None，则所有结构将重叠显示。默认为 None。
        show_sidechains (bool, optional): 是否显示侧链。默认为 False。
        show_mainchains (bool, optional): 是否显示主干。默认为 False。
        color_scheme (Literal["plddt", "rainbow", "chain"], optional): 
            颜色方案。默认为 "plddt"。

    Returns:
        py3Dmol.view: 可显示的可视化视图对象。
        
    Raises:
        FileNotFoundError: 如果任何一个 pdb_path 不存在。
        ValueError: 如果网格大小不足以容纳所有模型。
    """
    if grid_layout:
        rows, cols = grid_layout
        if len(pdb_paths) > rows * cols:
            raise ValueError(f"网格大小 ({rows}x{cols}) 不足以容纳 {len(pdb_paths)} 个模型。")
        view = py3Dmol.view(width=cols*400, height=rows*400, viewergrid=grid_layout, js='https://3dmol.org/build/3Dmol.js')
    else:
        # 如果没有提供网格布局，则退化为重叠显示
        view = py3Dmol.view(width=800, height=600, js='https://3dmol.org/build/3Dmol.js')

    # 循环加载和设置每个模型
    for i, pdb_path in enumerate(pdb_paths):
        pdb_file = Path(pdb_path)
        if not pdb_file.exists():
            raise FileNotFoundError(f"错误：找不到文件 '{pdb_path}'")
        
        # 计算当前模型在网格中的位置
        viewer_pos = (i // grid_layout[1], i % grid_layout[1]) if grid_layout else None
        
        # 将模型添加到指定的网格单元中
        view.addModel(pdb_file.read_text(), pdb_file.suffix[1:], viewer=viewer_pos)

        # 为当前网格单元中的模型设置样式
        if color_scheme == "plddt":
            view.setStyle({'cartoon': {'colorscheme': {'prop': 'b', 'gradient': 'roygb', 'min': 50, 'max': 90}}}, viewer=viewer_pos)
        elif color_scheme == "rainbow":
            view.setStyle({'cartoon': {'color': 'spectrum'}}, viewer=viewer_pos)
        elif color_scheme == "chain":
            view.setStyle({'cartoon': {'colorscheme': 'chain'}}, viewer=viewer_pos)

        if show_sidechains:
            backbone_atoms = ['C', 'O', 'N']
            view.addStyle({'and': [{'resn': ["GLY", "PRO"], 'invert': True}, {'atom': backbone_atoms, 'invert': True}]},
                          {'stick': {'colorscheme': 'lightgreyCarbon', 'radius': 0.3}}, viewer=viewer_pos)
            view.addStyle({'and': [{'resn': "GLY"}, {'atom': 'CA'}]},
                          {'sphere': {'colorscheme': 'lightgreyCarbon', 'radius': 0.3}}, viewer=viewer_pos)
            view.addStyle({'and': [{'resn': "PRO"}, {'atom': ['C', 'O'], 'invert': True}]},
                          {'stick': {'colorscheme': 'lightgreyCarbon', 'radius': 0.3}}, viewer=viewer_pos)
            
        if show_mainchains:
            backbone_atoms = ['C', 'O', 'N', 'CA']
            view.addStyle({'atom': backbone_atoms}, {'stick': {'colorscheme': 'lightgreyCarbon', 'radius': 0.3}}, viewer=viewer_pos)

        # 缩放以适应当前网格单元中的模型
        view.zoomTo(viewer=viewer_pos)

    return view