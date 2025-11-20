#!/usr/bin/env python3
"""
ldg_to_vtu.py
从 ldg.txt 解析 region 指令，生成每个区域的长方体 VTU 文件，便于可视化。
用法：python3 ldg_to_vtu.py finfet/ldg.txt output.vtu
"""

import sys
from pathlib import Path

def parse_regions(ldg_path):
    regions = []
    with open(ldg_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if parts[0].lower() == "region" and len(parts) >= 8:
                x1, x2 = float(parts[1]), float(parts[2])
                y1, y2 = float(parts[3]), float(parts[4])
                z1, z2 = float(parts[5]), float(parts[6])
                mat = parts[7]
                regions.append((x1, x2, y1, y2, z1, z2, mat))
    return regions

def write_vtu(regions, out_path):
    # 简单起见，每个 region 生成一个 hexahedron，点不去重
    points = []
    cells = []
    materials = []
    material_ids = []
    mat_to_id = {}

    for ridx, (x1, x2, y1, y2, z1, z2, mat) in enumerate(regions):
        pts = [
            (x1, y1, z1), (x2, y1, z1), (x2, y2, z1), (x1, y2, z1),  # 底面
            (x1, y1, z2), (x2, y1, z2), (x2, y2, z2), (x1, y2, z2),  # 顶面
        ]
        base = len(points)
        points.extend(pts)
        # VTK_HEXAHEDRON = cell type 12
        cells.append([base + i for i in range(8)])
        materials.append(mat)
        if mat not in mat_to_id:
            mat_to_id[mat] = len(mat_to_id)
        material_ids.append(mat_to_id[mat])

    with open(out_path, "w", encoding="utf-8") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        f.write('  <UnstructuredGrid>\n')
        f.write(f'    <Piece NumberOfPoints="{len(points)}" NumberOfCells="{len(cells)}">\n')
        # 点坐标
        f.write('      <Points>\n')
        f.write('        <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n')
        for p in points:
            f.write(f'          {p[0]} {p[1]} {p[2]}\n')
        f.write('        </DataArray>\n')
        f.write('      </Points>\n')
        # 单元连接
        f.write('      <Cells>\n')
        f.write('        <DataArray type="Int32" Name="connectivity" format="ascii">\n')
        for c in cells:
            f.write('          ' + ' '.join(map(str, c)) + '\n')
        f.write('        </DataArray>\n')
        f.write('        <DataArray type="Int32" Name="offsets" format="ascii">\n')
        offset = 0
        for c in cells:
            offset += len(c)
            f.write(f'          {offset}\n')
        f.write('        </DataArray>\n')
        f.write('        <DataArray type="UInt8" Name="types" format="ascii">\n')
        for _ in cells:
            f.write('          12\n')  # VTK_HEXAHEDRON
        f.write('        </DataArray>\n')
        f.write('      </Cells>\n')
        # 区域材料 ID 作为 CellData，方便配色
        f.write('      <CellData Scalars="material_id">\n')
        f.write('        <DataArray type="Int32" Name="material_id" format="ascii">\n')
        for mid in material_ids:
            f.write(f'          {mid}\n')
        f.write('        </DataArray>\n')
        f.write('      </CellData>\n')

        # 将材料名称表写入 FieldData 方便查阅
        f.write('      <FieldData>\n')
        f.write('        <DataArray type="String" Name="material_names" format="ascii">\n')
        # 按 ID 顺序写出名称
        for name, mid in sorted(mat_to_id.items(), key=lambda x: x[1]):
            f.write(f'          {name}\n')
        f.write('        </DataArray>\n')
        f.write('      </FieldData>\n')
        f.write('    </Piece>\n')
        f.write('  </UnstructuredGrid>\n')
        f.write('</VTKFile>\n')

def main():
    if len(sys.argv) != 3:
        print("用法: python3 ldg_to_vtu.py <ldg.txt> <output.vtu>")
        sys.exit(1)
    ldg_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2])
    regions = parse_regions(ldg_path)
    if not regions:
        print("没有找到 region 行，检查 ldg.txt")
        sys.exit(1)
    write_vtu(regions, out_path)
    print(f"写出 {len(regions)} 个区域到 {out_path}")

if __name__ == "__main__":
    main()
