import matplotlib.pyplot as plt
import matplotlib.patches as patches
from collections import defaultdict, deque
import re
import itertools  # 用于生成线段组合

def read_input_file(input_file):
    with open(input_file, 'r', encoding='utf-8', errors='ignore') as file:
        lines = file.readlines()
    
    dimX, dimY = 0, 0
    source = None
    sinks = []
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith('.dimx'):
            dimX = int(line.split()[1])
        elif line.startswith('.dimy'):
            dimY = int(line.split()[1])
        elif line.startswith('.e'):
            break
        elif line.startswith('.'):
            continue
        else:
            coords = list(map(float, line.strip().split()))
            if len(coords) >= 2:
                if source is None:
                    source = (coords[0], coords[1])  # 第一个坐标为源点
                else:
                    sinks.append((coords[0], coords[1]))  # 其余为汇点
            else:
                print(f"无效的坐标行: {line}")
    return dimX, dimY, source, sinks

def read_output_file(output_file):
    with open(output_file, 'r', encoding='utf-8', errors='ignore') as file:
        lines = file.readlines()

    segments = []
    num_segments = 0
    dimX, dimY = 0, 0

    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if line.startswith('.l'):
            num_segments = int(line.split()[1])
            print(f"總線段數量: {num_segments}")
        elif line.startswith('.dimx'):
            dimX = int(line.split()[1])
            print(f"晶片在 X 軸的尺寸: {dimX}")
        elif line.startswith('.dimy'):
            dimY = int(line.split()[1])
            print(f"晶片在 Y 軸的尺寸: {dimY}")
        elif line.startswith('.e'):
            break
        elif line.startswith('.'):
            continue
        else:
            coords = list(map(int, line.split()))
            if len(coords) == 4:
                x1, y1, x2, y2 = coords
                segments.append(((x1, y1), (x2, y2)))
            else:
                print(f"無效的線段行: {line}")

    return segments

def calculate_times(source_point, sink_points, segments):
    point_id = {}
    id_point = {}
    current_id = 0

    all_points = set()
    for segment in segments:
        all_points.add(segment[0])
        all_points.add(segment[1])
    if source_point:
        all_points.add(source_point)
    all_points.update(sink_points)

    for point in all_points:
        point_id[point] = current_id
        id_point[current_id] = point
        current_id += 1

    adj_list = defaultdict(list)
    for segment in segments:
        p1_id = point_id[segment[0]]
        p2_id = point_id[segment[1]]
        distance = abs(segment[0][0] - segment[1][0]) + abs(segment[0][1] - segment[1][1])
        adj_list[p1_id].append((p2_id, distance))
        adj_list[p2_id].append((p1_id, distance))

    if source_point is None:
        print("未提供源點，無法計算時間。")
        return 0, 0, sink_points  # 返回所有匯點為無法連接

    source_id = point_id[source_point]
    times = []
    unconnected_points = []

    for sink in sink_points:
        sink_id = point_id[sink]
        path_length = bfs(adj_list, source_id, sink_id)
        if path_length == float('inf'):
            print(f"警告：无法找到从源点 {source_point} 到汇点 {sink} 的路径。")
            times.append(0)
            unconnected_points.append(sink)
        else:
            times.append(path_length)

    valid_times = [t for t in times if t > 0]
    T_max = max(valid_times) if valid_times else 0
    T_min = min(valid_times) if valid_times else 0
    return T_max, T_min, unconnected_points

def bfs(adj_list, source_id, target_id):
    queue = deque([(source_id, 0)])  # (当前节点ID, 路径长度)
    visited = set()
    visited.add(source_id)

    while queue:
        current_id, total_length = queue.popleft()

        if current_id == target_id:
            return total_length

        for neighbor_id, length in adj_list[current_id]:
            if neighbor_id not in visited:
                visited.add(neighbor_id)
                queue.append((neighbor_id, total_length + length))

    return float('inf')

def calculate_wire_length(segments):
    total_length = 0
    for segment in segments:
        length = abs(segment[0][0] - segment[1][0]) + abs(segment[0][1] - segment[1][1])
        total_length += length
    return total_length

def visualize_cts_with_unconnected(dimX, dimY, source_point, sink_points, unconnected_points, segments, contained_pairs):
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.add_patch(patches.Rectangle((0, 0), dimX, dimY, linewidth=2, edgecolor='black', facecolor='none'))

    if source_point:
        ax.plot(source_point[0], source_point[1], 'ro', markersize=8, label='Source')

    if sink_points:
        sink_x = [p[0] for p in sink_points]
        sink_y = [p[1] for p in sink_points]
        ax.plot(sink_x, sink_y, 'bo', markersize=6, label='Sinks')

    if unconnected_points:
        unconnected_x = [p[0] for p in unconnected_points]
        unconnected_y = [p[1] for p in unconnected_points]
        ax.plot(unconnected_x, unconnected_y, 'yo', markersize=6, label='Unconnected Sinks')

    # 绘制所有线段
    for segment in segments:
        x_values = [segment[0][0], segment[1][0]]
        y_values = [segment[0][1], segment[1][1]]
        ax.plot(x_values, y_values, 'g-', linewidth=2)

    # 高亮显示完全包含的线段
    for ((i, seg_large), (j, seg_small)) in contained_pairs:
        # 绘制大的线段为红色虚线
        x_values_large = [seg_large[0][0], seg_large[1][0]]
        y_values_large = [seg_large[0][1], seg_large[1][1]]
        ax.plot(x_values_large, y_values_large, 'r--', linewidth=2, label='Containing Segment' if i == contained_pairs[0][0][0] else "")

        # 绘制小的线段为蓝色虚线
        x_values_small = [seg_small[0][0], seg_small[1][0]]
        y_values_small = [seg_small[0][1], seg_small[1][1]]
        ax.plot(x_values_small, y_values_small, 'b--', linewidth=2, label='Contained Segment' if i == contained_pairs[0][0][0] else "")

    title = f"Clock Tree Synthesis"
    ax.set_title(title, fontsize=14)
    ax.set_xlabel("X-axis", fontsize=12)
    ax.set_ylabel("Y-axis", fontsize=12)
    ax.legend(loc='upper right')

    ax.set_xlim(-10, dimX + 10)
    ax.set_ylim(-10, dimY + 20)

    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.show()

# 新增：检查线段完全包含的函数
def check_contained_segments(segments):
    contained_pairs = []
    for (i, seg_large), (j, seg_small) in itertools.combinations(enumerate(segments), 2):
        if segment_contains(seg_large, seg_small):
            contained_pairs.append(((i, seg_large), (j, seg_small)))
        elif segment_contains(seg_small, seg_large):
            contained_pairs.append(((j, seg_small), (i, seg_large)))
    return contained_pairs

def segment_contains(seg_large, seg_small):
    # 检查 seg_small 的两个端点是否都在 seg_large 上
    return on_segment(seg_large[0], seg_small[0], seg_large[1]) and on_segment(seg_large[0], seg_small[1], seg_large[1])

def on_segment(p, q, r):
    if min(p[0], r[0]) <= q[0] <= max(p[0], r[0]) and \
       min(p[1], r[1]) <= q[1] <= max(p[1], r[1]):
        return True
    return False

def main(input_file, output_file, W_FLUTE):
    dimX, dimY, source_input, sinks_input = read_input_file(input_file)
    segments = read_output_file(output_file)

    # 新增：检查线段完全包含的功能
    contained_pairs = check_contained_segments(segments)
    if contained_pairs:
        print(f"发现 {len(contained_pairs)} 对完全包含的线段：")
        for ((i, seg_large), (j, seg_small)) in contained_pairs:
            print(f"  线段 {i} {seg_large} 完全包含 线段 {j} {seg_small}")
    else:
        print("未发现完全包含的线段。")

    # 如果需要同时保留之前的重叠检查，可以添加相关代码
    # overlapping_pairs = check_overlapping_segments(segments)
    # if overlapping_pairs:
    #     print(f"发现 {len(overlapping_pairs)} 对重叠的线段：")
    #     for ((i, seg1), (j, seg2)) in overlapping_pairs:
    #         print(f"  线段 {i} {seg1} 与 线段 {j} {seg2} 重叠")
    # else:
    #     print("未发现重叠的线段。")

    T_max, T_min, unconnected_points = calculate_times(source_input, sinks_input, segments)
    W_cts = calculate_wire_length(segments)
    skew_ratio = T_max / T_min if T_min != 0 else float('inf')
    wire_length_ratio = W_cts / W_FLUTE if W_FLUTE != 0 else float('inf')

    print(f"T_max: {T_max}, T_min: {T_min}, Skew Ratio: {skew_ratio}")
    print(f"W_CTS: {W_cts}, W_FLUTE: {W_FLUTE}, Wire Length Ratio: {wire_length_ratio}")

    visualize_cts_with_unconnected(dimX, dimY, source_input, sinks_input, unconnected_points, segments, contained_pairs)

if __name__ == '__main__':
    input_file = 'case3.cts'  
    output_file = 'ouput03.cts'  # 您提供的檔案名稱
    W_FLUTE = 155  # 示例线长
    main(input_file, output_file, W_FLUTE)
