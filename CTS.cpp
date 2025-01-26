#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <cmath>
#include <climits>
#include <list>
#include <set>
#include <map>
#include <string>
#include <utility>

using namespace std;

// 定義 Point 結構
struct Point {
    int x, y;
    int id;  // 每個點的唯一識別碼

    // 用於比較兩個點是否相同
    bool operator==(const Point& other) const {
        return x == other.x && y == other.y;
    }

    // 新增建構函式，接受 x, y 和 id
    Point(int xVal, int yVal, int idVal) : x(xVal), y(yVal), id(idVal) {}

    // 預設建構函式
    Point() : x(0), y(0), id(-1) {}  // 可以提供預設值，避免未初始化
};

// 定義 Region 結構
struct Region {
    int minX, maxX, minY, maxY;
    vector<Point> sinks;
    Point massPoint;
    bool hasMassPoint;

    // 建構函式，接受5個引數
    Region(int min_x, int max_x, int min_y, int max_y, const Point& mass_point)
        : minX(min_x), maxX(max_x), minY(min_y), maxY(max_y),
          massPoint(mass_point), hasMassPoint(false) {}
};

// 定義 LargeRegion 結構
struct LargeRegion {
    int minX, maxX, minY, maxY;
    vector<Point> massPoints; // 包含的小區域質點
    Point massPoint;
    bool hasMassPoint;

    // 建構函式，接受5個引數
    LargeRegion(int min_x, int max_x, int min_y, int max_y, const Point& mass_point)
        : minX(min_x), maxX(max_x), minY(min_y), maxY(max_y),
          massPoint(mass_point), hasMassPoint(false) {}
};

// 自訂 hash 函式來支援 pair<int, int> 
struct hash_pair {
    template <class T1, class T2>
    size_t operator()(const pair<T1, T2>& p) const {
        auto hash1 = hash<T1>{}(p.first);
        auto hash2 = hash<T2>{}(p.second);
        return hash1 ^ (hash2 << 1); // XOR
    }
};

// 函式：讀取 CTS Input 檔案
vector<Point> readCTSInput(const string& filename, int& dimX, int& dimY, Point& source, vector<Point>& sinks) {
    vector<Point> points;
    ifstream inFile(filename);
    if (!inFile.is_open()) {
        cerr << "無法開啟輸入檔案: " << filename << endl;
        exit(1);
    }
    string line;
    int pointIndex = 0;

    while (getline(inFile, line)) {
        stringstream ss(line);
        if (line.empty()) continue;
        if (line[0] == '.') {
            string type;
            ss >> type;
            if (type == ".dimx") ss >> dimX;
            else if (type == ".dimy") ss >> dimY;
        }
        else if (line[0] != '.' && line[0] != 'e') {
            Point point;
            ss >> point.x >> point.y;
            point.id = pointIndex++;
            if (pointIndex == 1) {
                source = point;  // 第一個點作為源點
            }
            else {
                sinks.push_back(point);  // 其餘的點作為接收點
            }
            points.push_back(point);
        }
    }
    inFile.close();
    return points;
}

// 函式：計算質點(平均位置)
Point calculateMassPoint(const vector<Point>& points, int& nextId) {
    if (points.empty()) {
        return Point{ 0, 0, -1 };
    }
    double sumX = 0, sumY = 0;
    for (const auto& p : points) {
        sumX += p.x;
        sumY += p.y;
    }
    int massX = static_cast<int>(round(sumX / points.size()));
    int massY = static_cast<int>(round(sumY / points.size()));
    Point massPoint = Point{ massX, massY, nextId++ };
    return massPoint;
}

// 函式：將所有線段映射到2D陣列，並標記Sink點
vector<vector<int>> mapTo2DArray(const vector<pair<Point, Point>>& segments, int dimX, int dimY, const vector<Point>& sinks) {
    vector<vector<int>> grid(dimY + 1, vector<int>(dimX + 1, 0));

    // 映射線段到 grid
    for (const auto& seg : segments) {
        Point p1 = seg.first;
        Point p2 = seg.second;

        if (p1.x == p2.x) { // 垂直線
            int x = p1.x;
            int startY = min(p1.y, p2.y);
            int endY = max(p1.y, p2.y);
            for (int y = startY; y <= endY; ++y) {
                grid[y][x] = 1;
            }
        }
        else if (p1.y == p2.y) { // 水平線
            int y = p1.y;
            int startX = min(p1.x, p2.x);
            int endX = max(p1.x, p2.x);
            for (int x = startX; x <= endX; ++x) {
                grid[y][x] = 1;
            }
        }
        else { // 斜線不處理
   
        }
    }

    // 將 Sink 點標記在 grid 上
    for (const auto& sink : sinks) {
        if (sink.x >= 0 && sink.x <= dimX && sink.y >= 0 && sink.y <= dimY) {
            grid[sink.y][sink.x] = 1;
        }
    }

    return grid;
}

// 函式：檢查線段是否為水平或垂直線段
bool isHorizontal(const pair<Point, Point>& seg) {
    return seg.first.y == seg.second.y;
}

bool isVertical(const pair<Point, Point>& seg) {
    return seg.first.x == seg.second.x;
}

// 函式：檢測水平與垂直線段是否相交
bool checkIntersection(const pair<Point, Point>& horizontal, const pair<Point, Point>& vertical, Point& intersection) {
    if (isHorizontal(horizontal) && isVertical(vertical)) {
        if (vertical.first.x >= min(horizontal.first.x, horizontal.second.x) &&
            vertical.first.x <= max(horizontal.first.x, horizontal.second.x) &&
            horizontal.first.y >= min(vertical.first.y, vertical.second.y) &&
            horizontal.first.y <= max(vertical.first.y, vertical.second.y)) {
            intersection = Point(vertical.first.x, horizontal.first.y, -1);  // 使用完整初始化
            return true;
        }
    }
    return false;
}

// 函式：輸出自訂格式的線段資料
void outputCustomFormat(const string& filename, const vector<pair<Point, Point>>& segments, int dimX, int dimY) {
    ofstream outFile(filename);
    if (!outFile.is_open()) {
        cerr << "無法開啟輸出檔案: " << filename << endl;
        return;
    }

    // 輸出線段數量和晶片尺寸
    outFile << ".l " << segments.size() << "\n";
    outFile << ".dimx " << dimX << "\n";
    outFile << ".dimy " << dimY << "\n";

    // 輸出每個線段的座標
    for (const auto& seg : segments) {
        outFile << seg.first.x << " " << seg.first.y << " " << seg.second.x << " " << seg.second.y << "\n";
    }

    // 結尾標記
    outFile << ".e\n";

    outFile.close();
}

// 主輸出函式，處理分段和輸出
void outputGnuplot(const string& filename, vector<pair<Point, Point>>& segments,
    const Point& source, const vector<Point>& sinks, int dimX, int dimY,
    const vector<Point>& failedPoints) {

    // 使用 unordered_map 來記錄每個點的交點
    unordered_map<pair<int, int>, vector<Point>, hash_pair> intersectionMap;

    // 遍歷每個線段，檢查是否有交點，並記錄這些交點
    for (const auto& seg1 : segments) {
        if (isHorizontal(seg1)) {
            // 找垂直交點
            for (const auto& seg2 : segments) {
                if (isVertical(seg2)) {
                    Point intersection;
                    if (checkIntersection(seg1, seg2, intersection)) {
                        // 記錄交點
                        intersectionMap[{seg1.first.x, seg1.first.y}].push_back(intersection);
                        intersectionMap[{seg2.first.x, seg2.first.y}].push_back(intersection);
                    }
                }
            }
        }
    }

    // 處理與接收點有關的拆分
    for (const auto& sink : sinks) {
        for (auto& seg : segments) {
            if ((isHorizontal(seg) && sink.y == seg.first.y &&
                sink.x > min(seg.first.x, seg.second.x) && sink.x < max(seg.first.x, seg.second.x)) ||
                (isVertical(seg) && sink.x == seg.first.x &&
                    sink.y > min(seg.first.y, seg.second.y) && sink.y < max(seg.first.y, seg.second.y))) {
                // 在線段中插入接收點
                intersectionMap[{seg.first.x, seg.first.y}].push_back(sink);
            }
        }
    }

    // 現在進行線段的拆分
    vector<pair<Point, Point>> splitSegments;
    for (const auto& seg : segments) {
        Point p1 = seg.first;
        Point p2 = seg.second;

        vector<Point> splitPoints = { p1 };  // 初始化分段點數組從 p1 開始

        // 檢查該線段的分段點
        auto key = make_pair(p1.x, p1.y);
        if (intersectionMap.find(key) != intersectionMap.end()) {
            const auto& intersections = intersectionMap[key];
            for (const auto& interPoint : intersections) {
                // 如果這個交點在線段範圍內加入分段
                if (isHorizontal(seg) && interPoint.x > min(p1.x, p2.x) && interPoint.x < max(p1.x, p2.x)) {
                    splitPoints.push_back(interPoint);
                }
                else if (isVertical(seg) && interPoint.y > min(p1.y, p2.y) && interPoint.y < max(p1.y, p2.y)) {
                    splitPoints.push_back(interPoint);
                }
            }
        }

        splitPoints.push_back(p2);  // 最後加入終點 p2

        // 將分段點組成線段
        sort(splitPoints.begin(), splitPoints.end(), [](const Point& a, const Point& b) {
            return (a.x == b.x) ? a.y < b.y : a.x < b.x;
            });

        // 將分段結果放入splitSegments
        for (size_t i = 0; i < splitPoints.size() - 1; ++i) {
            splitSegments.emplace_back(splitPoints[i], splitPoints[i + 1]);
        }
    }

    // 輸出
    outputCustomFormat(filename, splitSegments, dimX, dimY);
}

//處理分段和輸出
void outputResults(const string& ctsOutputFilename,
    const vector<pair<Point, Point>>& segments, const Point& source, const vector<Point>& sinks,
    const vector<Point>& largeRegionMassPoints, const vector<Point>& regionMassPoints,
    int dimX, int dimY) {

    // 輸出 CTS 輸出檔案
    outputGnuplot(ctsOutputFilename, const_cast<vector<pair<Point, Point>>&>(segments), source, sinks, dimX, dimY, vector<Point>()); // failedPoints 未使用
    //outputCustomFormat(ctsOutputFilename, segments, dimX, dimY);
    cout << "CTS 輸出已寫入至 '" << ctsOutputFilename << "'。\n";

}

// 主程式
int main(int argc, char* argv[]) {
    // 檢查命令行參數
    if (argc != 3) {
        cerr << "使用方式: " << argv[0] << " <input_file> <output_file>" << endl;
        return 1;
    }

    string inputFile = argv[1];  // 第一個引數是輸入檔案
    string outputFile = argv[2]; // 第二個引數是輸出檔案

    int dimX = 0, dimY = 0;
    Point source;
    vector<Point> sinks;

    // 步驟 1: 讀取 CTS input 檔案，確保檔案存在且格式正確
    vector<Point> points = readCTSInput(inputFile, dimX, dimY, source, sinks);
    if (points.empty()) {
        cerr << "讀取輸入檔案時發生錯誤。" << endl;
        return 1;
    }

    // 初始化 nextId 為現有點的數量，以便分配給新點
    int nextId = static_cast<int>(points.size());

    // 步驟 2: 計算整個座標系統的總質點
    Point overallMassPoint = calculateMassPoint(points, nextId);

    // 步驟 3: 定義4個大區域 (2x2網格)
    int largeRegionCols = 2;
    int largeRegionRows = 2;
    int largeRegionWidth = (dimX + largeRegionCols - 1) / largeRegionCols; // 
    int largeRegionHeight = (dimY + largeRegionRows - 1) / largeRegionRows;

    vector<LargeRegion> largeRegions;
    for (int i = 0; i < largeRegionCols; ++i) {
        for (int j = 0; j < largeRegionRows; ++j) {
            int minX = i * largeRegionWidth;
            int maxX = min((i + 1) * largeRegionWidth - 1, dimX); // 保證邊界
            int minY = j * largeRegionHeight;
            int maxY = min((j + 1) * largeRegionHeight - 1, dimY); // 保證邊界

            // 計算大區域的中心點
            int centerX = (minX + maxX) / 2;
            int centerY = (minY + maxY) / 2;

            // 新增這個大區域
            largeRegions.emplace_back(minX, maxX, minY, maxY, Point(centerX, centerY, nextId++));
        }
    }

    // 步驟 4: 分配大區域的質點到各自的大區域並計算大區域質點
    vector<Point> largeRegionMassPoints;
    for (auto& largeRegion : largeRegions) {
        // 計算大區域的質點
        // 在這裡可以根據需要分配點到大區域，這裡暫時假設每個大區域的質點已經定義
        largeRegion.hasMassPoint = true;
        largeRegionMassPoints.push_back(largeRegion.massPoint);
    }

    vector<pair<Point, Point>> sourceSegs;
    vector<Point> failedPoints;
    vector<pair<Point, Point>> segments;
    
    // 步驟 5: 連接大區域質點到總質點，先垂直再水平，並加入重疊碰撞檢查
    // 這裡假設從每個大區域質點直接連接到總質點

    // 初始化 grid 並標記初始的 segments 和 sinks
    vector<vector<int>> grid = mapTo2DArray(segments, dimX, dimY, sinks);

    for (const auto& largeRegion : largeRegions) {
        Point currentPoint = largeRegion.massPoint;
        bool collision = false;  // 初始化碰撞標記
        vector<pair<Point, Point>> newSegs;

        // 先進行垂直連接：從大區域質點連接到總質點的 y 坐標
        int stepY = (overallMassPoint.y >= largeRegion.massPoint.y) ? 1 : -1;
        for (int y = largeRegion.massPoint.y; y != overallMassPoint.y; y += stepY) {
            int nextY = y + stepY;

            // 確認是否發生碰撞（從 grid 判斷碰撞）
            if (grid[nextY][largeRegion.massPoint.x] == 1) {
                collision = true;
                Point failedPoint = Point{ largeRegion.massPoint.x, nextY, -1 };
                newSegs.emplace_back(currentPoint, failedPoint);  // 保留碰撞前的線段
                currentPoint = failedPoint;
                break;  // 停止垂直連接
            }

            // 沒有碰撞，則繼續連接
            Point nextPoint = Point{ largeRegion.massPoint.x, nextY, -1 };
            newSegs.emplace_back(currentPoint, nextPoint);
            currentPoint = nextPoint;

            // 更新 grid
            grid[nextPoint.y][nextPoint.x] = 1;
        }

        // 水平連接：從垂直終點連接到總質點的 x 坐標
        if (!collision && currentPoint.y == overallMassPoint.y && currentPoint.x != overallMassPoint.x) {
            int stepX = (overallMassPoint.x >= currentPoint.x) ? 1 : -1;
            for (int x = currentPoint.x; x != overallMassPoint.x; x += stepX) {
                int nextX = x + stepX;

                // 確認是否發生碰撞（從 grid 判斷碰撞）
                if (grid[currentPoint.y][nextX] == 1) {
                    collision = true;
                    Point failedPoint = Point{ nextX, currentPoint.y, -1 };
                    newSegs.emplace_back(currentPoint, failedPoint);  // 保留碰撞前的線段
                    currentPoint = failedPoint;
                    break;  // 停止水平連接
                }

                // 沒有碰撞，則繼續連接
                Point nextPoint = Point{ nextX, currentPoint.y, -1 };
                newSegs.emplace_back(currentPoint, nextPoint);
                currentPoint = nextPoint;

                // 更新 grid
                grid[nextPoint.y][nextPoint.x] = 1;
            }
        }

        // 如果發生了碰撞，可以嘗試替代路徑或記錄失敗點
        if (collision) {
            failedPoints.push_back(currentPoint);  // 記錄碰撞的失敗點
        }

        // 合併相鄰的線段
        vector<pair<Point, Point>> mergedSegs;
        if (!newSegs.empty()) {
            Point start = newSegs[0].first;
            Point end = newSegs[0].second;

            for (size_t i = 1; i < newSegs.size(); ++i) {
                Point nextStart = newSegs[i].first;
                Point nextEnd = newSegs[i].second;

                // 檢查是否可以合併
                if ((start.x == end.x && end.x == nextStart.x && nextStart.x == nextEnd.x) ||
                    (start.y == end.y && end.y == nextStart.y && nextStart.y == nextEnd.y)) {
                    end = nextEnd;
                }
                else {
                    mergedSegs.emplace_back(start, end);
                    start = nextStart;
                    end = nextEnd;
                }
            }
            mergedSegs.emplace_back(start, end);
        }

        // 插入合併後的線段
        segments.insert(segments.end(), mergedSegs.begin(), mergedSegs.end());
        mergedSegs.clear();
        newSegs.clear();
    }

    // 步驟 6: 定義16個小區域 (4x4網格)
    int smallRegionCols = 4;
    int smallRegionRows = 4;
    int smallRegionWidth = (dimX + smallRegionCols - 1) / smallRegionCols; // 確保除不盡時，區域不會少
    int smallRegionHeight = (dimY + smallRegionRows - 1) / smallRegionRows;

    vector<Region> regions;
    for (int i = 0; i < smallRegionCols; ++i) {
        for (int j = 0; j < smallRegionRows; ++j) {
            // 定義區域的最小和最大 X, Y
            int minX = i * smallRegionWidth;
            int maxX = min((i + 1) * smallRegionWidth - 1, dimX); // 保證區域邊界不超過系統範圍
            int minY = j * smallRegionHeight;
            int maxY = min((j + 1) * smallRegionHeight - 1, dimY); // 同樣保證區域邊界的範圍

            // 計算區域的中心點 (以該區域的質點作為中心)
            int centerX = (minX + maxX) / 2;
            int centerY = (minY + maxY) / 2;

            // 新增這個區域
            regions.emplace_back(minX, maxX, minY, maxY, Point(centerX, centerY, nextId++));
        }
    }

    // 步驟 7: 分配 sinks 到各個小區域並計算小區域質點
    vector<Point> regionMassPoints;
    for (auto& region : regions) {
        // 分配 sinks 到區域
        for (const auto& sink : sinks) {
            if (sink.x >= region.minX && sink.x <= region.maxX &&
                sink.y >= region.minY && sink.y <= region.maxY) {
                region.sinks.push_back(sink);
            }
        }

        // 計算質點如果有 sinks
        if (!region.sinks.empty()) {
            region.massPoint = calculateMassPoint(region.sinks, nextId);
            region.hasMassPoint = true;
            regionMassPoints.push_back(region.massPoint);
        }
    }

    // 步驟 8: 連接大區域的質點與小區域的質點，先水平再垂直，從小區域質點往大區域質點方向進行碰撞檢查與合併相鄰線段
    for (const auto& largeRegion : largeRegions) {
        for (const auto& region : regions) {
            // 檢查小區域是否在大區域內
            if (region.minX >= largeRegion.minX && region.maxX <= largeRegion.maxX &&
                region.minY >= largeRegion.minY && region.maxY <= largeRegion.maxY) {
                if (region.hasMassPoint) {
                    // 初始化當前點為小區域質點
                    Point currentPoint = region.massPoint;
                    bool collision = false;  // 初始化碰撞標記
                    vector<pair<Point, Point>> regionSegs;

                    // 先水平連接：從小區域質點往大區域質點方向
                    int stepX = (largeRegion.massPoint.x >= region.massPoint.x) ? 1 : -1;
                    for (int x = region.massPoint.x; x != largeRegion.massPoint.x; x += stepX) {
                        int nextX = x + stepX;

                        // 確認是否發生碰撞（從 grid 判斷碰撞）
                        if (grid[region.massPoint.y][nextX] == 1) {
                            collision = true;
                            Point failedPoint = Point{ nextX, region.massPoint.y, -1 };
                            failedPoints.push_back(failedPoint);
                            regionSegs.emplace_back(currentPoint, failedPoint);
                            currentPoint = failedPoint;
                            break;  // 停止水平連接
                        }

                        // 沒有碰撞，則繼續連接
                        Point nextPoint = Point{ nextX, region.massPoint.y, -1 };
                        regionSegs.emplace_back(currentPoint, nextPoint);
                        currentPoint = nextPoint;

                        // 更新 grid
                        grid[nextPoint.y][nextPoint.x] = 1;
                    }

                    // 再垂直連接：從小區域質點往大區域質點方向
                    if (!collision && currentPoint.x == largeRegion.massPoint.x && currentPoint.y != largeRegion.massPoint.y) {
                        int stepY = (largeRegion.massPoint.y >= currentPoint.y) ? 1 : -1;
                        for (int y = currentPoint.y; y != largeRegion.massPoint.y; y += stepY) {
                            int nextY = y + stepY;

                            // 確認是否發生碰撞（從 grid 判斷碰撞）
                            if (grid[nextY][currentPoint.x] == 1) {
                                collision = true;
                                Point failedPoint = Point{ currentPoint.x, nextY, -1 };
                                failedPoints.push_back(failedPoint);
                                regionSegs.emplace_back(currentPoint, failedPoint);
                                currentPoint = failedPoint;
                                break;  // 停止垂直連接
                            }

                            // 沒有碰撞，則繼續連接
                            Point nextPoint = Point{ currentPoint.x, nextY, -1 };
                            regionSegs.emplace_back(currentPoint, nextPoint);
                            currentPoint = nextPoint;

                            // 更新 grid
                            grid[nextPoint.y][nextPoint.x] = 1;
                        }
                    }

                    // 如果發生了碰撞，可以嘗試替代路徑或記錄失敗點
                    if (collision) {
                        failedPoints.push_back(currentPoint);  // 記錄碰撞的失敗點
                    }

                    // 合併相鄰的線段
                    vector<pair<Point, Point>> mergedSegs;
                    if (!regionSegs.empty()) {
                        Point start = regionSegs[0].first;
                        Point end = regionSegs[0].second;

                        for (size_t i = 1; i < regionSegs.size(); ++i) {
                            Point nextStart = regionSegs[i].first;
                            Point nextEnd = regionSegs[i].second;

                            // 檢查是否可以合併
                            if ((start.x == end.x && end.x == nextStart.x && nextStart.x == nextEnd.x) ||
                                (start.y == end.y && end.y == nextStart.y && nextStart.y == nextEnd.y)) {
                                end = nextEnd;
                            }
                            else {
                                mergedSegs.emplace_back(start, end);
                                start = nextStart;
                                end = nextEnd;
                            }
                        }
                        mergedSegs.emplace_back(start, end);
                    }

                    // 插入合併後的線段
                    segments.insert(segments.end(), mergedSegs.begin(), mergedSegs.end());
                    grid = mapTo2DArray(segments, dimX, dimY, sinks); // 更新 grid
                    mergedSegs.clear();
                    regionSegs.clear();
                }
            }
        }
    }

    // 步驟 9: 連接小區域質點與 sinks，先垂直再水平，從接收點(sink)往小區域質點方向進行碰撞檢查與合併相鄰線段
    for (const auto& region : regions) {
        if (region.hasMassPoint) {
            for (const auto& sink : region.sinks) {
                Point currentPoint = sink;  // 初始化當前點為sink
                bool collision = false;  // 初始化碰撞標記
                vector<pair<Point, Point>> sinkSegs;

                // 先垂直連接：從接收點往小區域質點方向
                int stepY = (region.massPoint.y >= sink.y) ? 1 : -1;
                for (int y = sink.y; y != region.massPoint.y; y += stepY) {
                    int nextY = y + stepY;
                    if (grid[nextY][sink.x] == 1) {
                        collision = true;
                        Point failedPoint = Point{ sink.x, nextY, -1 };
                        failedPoints.push_back(failedPoint);
                        sinkSegs.emplace_back(currentPoint, failedPoint);
                        currentPoint = failedPoint;
                        break;  // 停止垂直連接
                    }
                    Point nextPoint = Point{ sink.x, nextY, -1 };
                    sinkSegs.emplace_back(currentPoint, nextPoint);
                    currentPoint = nextPoint;

                    // 更新 grid
                    grid[nextPoint.y][nextPoint.x] = 1;
                }

                // 再水平連接：從接收點往小區域質點方向
                if (!collision && currentPoint.y == region.massPoint.y && currentPoint.x != region.massPoint.x) {
                    int stepX = (region.massPoint.x >= currentPoint.x) ? 1 : -1;
                    for (int x = currentPoint.x; x != region.massPoint.x; x += stepX) {
                        int nextX = x + stepX;
                        if (grid[currentPoint.y][nextX] == 1) {
                            collision = true;
                            Point failedPoint = Point{ nextX, currentPoint.y, -1 };
                            failedPoints.push_back(failedPoint);
                            sinkSegs.emplace_back(currentPoint, failedPoint);
                            currentPoint = failedPoint;
                            break;  // 停止水平連接
                        }
                        Point nextPoint = Point{ nextX, currentPoint.y, -1 };
                        sinkSegs.emplace_back(currentPoint, nextPoint);
                        currentPoint = nextPoint;

                        // 更新 grid
                        grid[nextPoint.y][nextPoint.x] = 1;
                    }
                }

                // 如果發生了碰撞，可以嘗試替代路徑或記錄失敗點
                if (collision) {
                    failedPoints.push_back(currentPoint);  // 記錄碰撞的失敗點
                }

                // 合併相鄰的線段
                vector<pair<Point, Point>> mergedSegs;
                if (!sinkSegs.empty()) {
                    Point start = sinkSegs[0].first;
                    Point end = sinkSegs[0].second;

                    for (size_t i = 1; i < sinkSegs.size(); ++i) {
                        Point nextStart = sinkSegs[i].first;
                        Point nextEnd = sinkSegs[i].second;

                        // 檢查是否可以合併
                        if ((start.x == end.x && end.x == nextStart.x && nextStart.x == nextEnd.x) ||
                            (start.y == end.y && end.y == nextStart.y && nextStart.y == nextEnd.y)) {
                            end = nextEnd;
                        }
                        else {
                            mergedSegs.emplace_back(start, end);
                            start = nextStart;
                            end = nextEnd;
                        }
                    }
                    mergedSegs.emplace_back(start, end);
                }

                // 插入合併後的線段
                segments.insert(segments.end(), mergedSegs.begin(), mergedSegs.end());
                grid = mapTo2DArray(segments, dimX, dimY, sinks); // 更新 grid
                mergedSegs.clear();
                sinkSegs.clear();
            }
        }
    }

    // 步驟 10: 連接源點與總質點，並準確判斷碰撞點
    // 初始化當前點為源點
    Point currentPoint = source;

    // 嘗試一次性連接源點和總質點
    bool collision = false;

    // 檢查水平線上的碰撞
    int stepX = (overallMassPoint.x >= source.x) ? 1 : -1;
    for (int x = source.x; x != overallMassPoint.x; x += stepX) {
        if (grid[source.y][x] == 1) {
            collision = true;
            break;
        }
    }

    // 檢查垂直線上的碰撞
    if (!collision) {
        int stepY = (overallMassPoint.y >= source.y) ? 1 : -1;
        for (int y = source.y; y != overallMassPoint.y; y += stepY) {
            if (grid[y][overallMassPoint.x] == 1) {
                collision = true;
                break;
            }
        }
    }
    // 如果沒有碰撞，直接連接源點到總質點
    if (!collision) {
        sourceSegs.emplace_back(source, overallMassPoint);
        // 更新 grid
        grid[overallMassPoint.y][overallMassPoint.x] = 1;
    }
    else {
        // 水平連接
        for (int x = source.x; x != overallMassPoint.x; x += stepX) {
            int nextX = x + stepX;

            // 確認是否發生碰撞（從 grid 判斷碰撞）
            if (grid[source.y][nextX] == 1) {
                // 碰撞發生，記錄失敗點並停止連接
                Point failedPoint = Point{ nextX, source.y, -1 };
                failedPoints.push_back(failedPoint);
                sourceSegs.emplace_back(currentPoint, failedPoint); // 保留碰撞前的線段
                currentPoint = failedPoint;
                break; // 停止連接
            }

            // 沒有碰撞，則繼續連接
            Point nextPoint = Point{ nextX, source.y, -1 };
            sourceSegs.emplace_back(currentPoint, nextPoint);
            currentPoint = nextPoint;

            // 更新 grid
            grid[nextPoint.y][nextPoint.x] = 1;
        }

        // 垂直連接
        if (currentPoint.x == overallMassPoint.x && currentPoint.y != overallMassPoint.y) {
            int stepY = (overallMassPoint.y >= currentPoint.y) ? 1 : -1;
            for (int y = currentPoint.y; y != overallMassPoint.y; y += stepY) {
                int nextY = y + stepY;

                // 確認是否發生碰撞（從 grid 判斷碰撞）
                if (grid[nextY][currentPoint.x] == 1) {
                    // 碰撞發生，記錄失敗點並停止連接
                    Point failedPoint = Point{ currentPoint.x, nextY, -1 };
                    failedPoints.push_back(failedPoint);
                    sourceSegs.emplace_back(currentPoint, failedPoint); // 保留碰撞前的線段
                    currentPoint = failedPoint;
                    break; // 停止連接
                }

                // 沒有碰撞，則繼續連接
                Point nextPoint = Point{ currentPoint.x, nextY, -1 };
                sourceSegs.emplace_back(currentPoint, nextPoint);
                currentPoint = nextPoint;

                // 更新 grid
                grid[nextPoint.y][nextPoint.x] = 1;
            }
        }
    }
    // 合併相鄰的線段
    vector<pair<Point, Point>> mergedSegs;

    if (!sourceSegs.empty()) {
        Point start = sourceSegs[0].first;
        Point end = sourceSegs[0].second;

        for (size_t i = 1; i < sourceSegs.size(); ++i) {
            Point nextStart = sourceSegs[i].first;
            Point nextEnd = sourceSegs[i].second;

            // 檢查是否可以合併（相鄰且在同一直線上）
            if ((start.x == end.x && end.x == nextStart.x && nextStart.x == nextEnd.x) ||  // 垂直線
                (start.y == end.y && end.y == nextStart.y && nextStart.y == nextEnd.y)) {  // 水平線
                // 可以合併，更新 end 點
                end = nextEnd;
            }
            else {
                // 無法合併，先將當前線段加入 mergedSegs
                mergedSegs.emplace_back(start, end);
                // 更新 start 和 end 為下一條線段
                start = nextStart;
                end = nextEnd;
            }
        }

        mergedSegs.emplace_back(start, end);
    }

    // 插入合併後的線段到 segments
    segments.insert(segments.end(), mergedSegs.begin(), mergedSegs.end());
    grid = mapTo2DArray(segments, dimX, dimY, sinks); // 更新 grid
    mergedSegs.clear();

    // 收集所有質點
    vector<Point> allMassPoints;
    // 加入總質點
    allMassPoints.push_back(overallMassPoint);
    // 加入大區域質點
    allMassPoints.insert(allMassPoints.end(), largeRegionMassPoints.begin(), largeRegionMassPoints.end());
    // 加入小區域質點
    allMassPoints.insert(allMassPoints.end(), regionMassPoints.begin(), regionMassPoints.end());

    // 步驟 11: 合併生成並輸出以及將結果寫入 CTS 輸出檔案
    outputResults(outputFile, segments, source, sinks, largeRegionMassPoints, regionMassPoints, dimX, dimY);

    return 0;
}
