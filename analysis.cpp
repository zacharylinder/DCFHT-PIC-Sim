#include <windows.h>
#include <windowsx.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <climits>
#include <gdiplus.h>
#include <cmath>
#include <random>
#include <iomanip>
#include <cfloat>

#pragma comment (lib,"Gdiplus.lib")

using namespace Gdiplus;

// -----------------------------------------------------------------------------
// Global Constants and State
// -----------------------------------------------------------------------------
constexpr int SCROLL_BAR_WIDTH = 200;
constexpr int SCROLL_BAR_HEIGHT = 20;
constexpr int THUMB_WIDTH = 20;
constexpr int MARGIN_BOTTOM = 10;
constexpr int BUTTON_WIDTH = 100;
constexpr int BUTTON_BUFFER = 10;

ULONG_PTR g_gdiplusToken;
int dataset = 1;
bool drawParticlesFlag = true;
bool drawCellsFlag = false;
int scrollValue = 0; // now represents an index (0 to totalDatasets-1)
bool isDragging = false;
int dragStartX = 0;
int thumbStartPos = 0;
int computedScrollBarX = 0;
int computedScrollBarY = 0;
int optionsCells = 1;
bool showGridPoints = false;
int radio = 2;

// New global variable to store the total number of datasets.
int totalDatasets = 1;  // default to 1 if file reading fails

// -----------------------------------------------------------------------------
// Particle structure
// -----------------------------------------------------------------------------
struct Particle {
    int particle, step;
    double x, y, z;
};

// -----------------------------------------------------------------------------
// RAII wrapper for GDI objects to auto-delete
// -----------------------------------------------------------------------------
class GDIObject {
public:
    HGDIOBJ handle;
    explicit GDIObject(HGDIOBJ h) : handle(h) {}
    ~GDIObject() {
        if (handle)
            DeleteObject(handle);
    }
};

// -----------------------------------------------------------------------------
// GDI+ Initialization and Shutdown
// -----------------------------------------------------------------------------
void InitGDIPlus() {
    GdiplusStartupInput input;
    GdiplusStartup(&g_gdiplusToken, &input, nullptr);
}

void ShutdownGDIPlus() {
    GdiplusShutdown(g_gdiplusToken);
}

int countDatasets(const std::string& filename) {
    std::ifstream file(filename);
    if (!file)
        return 0;

    int count = 1;  // start at 1 for the first dataset
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty())
            ++count;
    }
    return count - 1;
}

// -----------------------------------------------------------------------------
// Draw Scrollbar: draws the track and thumb, and updates the dataset value.
// -----------------------------------------------------------------------------
void drawScrollBar(HDC hdc, int x, int y) {
    // Draw track
    {
        HBRUSH hBrush = CreateSolidBrush(RGB(180, 180, 180));
        RECT track = { x, y, x + SCROLL_BAR_WIDTH, y + SCROLL_BAR_HEIGHT };
        FillRect(hdc, &track, hBrush);
        DeleteObject(hBrush);
    }
    // Draw thumb using the calculated totalDatasets.
    // Ensure we don't divide by zero (when totalDatasets == 1, thumb stays at x)
    int thumbPos = x;
    if (totalDatasets > 1)
        thumbPos = x + (scrollValue * (SCROLL_BAR_WIDTH - THUMB_WIDTH)) / (totalDatasets - 1);
    // Draw thumb
    {
        HBRUSH hBrush = CreateSolidBrush(RGB(100, 100, 100));
        RECT thumb = { thumbPos, y, thumbPos + THUMB_WIDTH, y + SCROLL_BAR_HEIGHT };
        FillRect(hdc, &thumb, hBrush);
        DeleteObject(hBrush);
    }
    // Update dataset; we add 1 because scrollValue is 0-indexed.
    dataset = std::min(scrollValue + 1, totalDatasets);
}

// -----------------------------------------------------------------------------
// Draw Particles: reads "log.txt", computes the center, and draws each particle
// with a cached brush (color) per simulation step.
// -----------------------------------------------------------------------------
void drawParticles(HWND hwnd, HDC hdc, int viewCenterX, int viewCenterY) {
    std::ifstream logFile("log.txt");
    if (!logFile) {
        MessageBox(hwnd, "Failed to open log.txt", "Error", MB_OK | MB_ICONERROR);
        return;
    }

    std::vector<Particle> particles;
    Particle p;
    // Also track simulation bounds.
    double simMinX = DBL_MAX, simMaxX = -DBL_MAX;
    double simMinY = DBL_MAX, simMaxY = -DBL_MAX;
    while (logFile >> p.particle >> p.step >> p.x >> p.y >> p.z) {
        particles.push_back(p);
        simMinX = std::min(simMinX, p.x);
        simMaxX = std::max(simMaxX, p.x);
        simMinY = std::min(simMinY, p.y);
        simMaxY = std::max(simMaxY, p.y);
    }
    logFile.close();

    if (particles.empty())
        return;

    int count = static_cast<int>(particles.size());
    double xCenter = (simMinX + simMaxX) / 2;
    double yCenter = (simMinY + simMaxY) / 2;

    // Get client area dimensions to compute a dynamic scaling factor.
    RECT clientRect;
    GetClientRect(hwnd, &clientRect);
    int availWidth = clientRect.right - clientRect.left;
    int availHeight = clientRect.bottom - clientRect.top;
    double simWidth = simMaxX - simMinX;
    double simHeight = simMaxY - simMinY;
    double scale = 10000;  // fallback value
    if (simWidth > 0 && simHeight > 0) {
        // Scale so that simulation fits in 90% of the view.
        scale = std::min(availWidth * 0.75 / simWidth, availHeight * 0.75 / simHeight);
    }

    if (radio == 1)
    {
        std::unordered_map<int, COLORREF> stepColors;
        std::unordered_map<int, HBRUSH> brushCache;
        HPEN nullPen = CreatePen(PS_NULL, 1, RGB(255, 255, 255));
        HPEN oldPen = (HPEN)SelectObject(hdc, nullPen);
        std::mt19937 rng(std::random_device{}());
        std::uniform_int_distribution<int> dist(0, 255);

        for (const auto& part : particles) {
            if (stepColors.find(part.step) == stepColors.end()) {
                COLORREF col = RGB(dist(rng), dist(rng), dist(rng));
                stepColors[part.step] = col;
                brushCache[part.step] = CreateSolidBrush(col);
            }
            SelectObject(hdc, brushCache[part.step]);
            // Compute screen coordinates using the dynamic scale.
            int px = static_cast<int>(viewCenterX + (xCenter - part.x) * scale);
            int py = static_cast<int>(viewCenterY + (yCenter - part.y) * scale);
            Ellipse(hdc, px - 1, py - 1, px + 1, py + 1);
        }

        // Clean up brushes and restore pen.
        for (auto& pair : brushCache) {
            DeleteObject(pair.second);
        }
        SelectObject(hdc, oldPen);
        DeleteObject(nullPen);
    }
    else if (radio == 2)
    {
        // Transform particles to screen coordinates and determine bounds.
        std::vector<std::pair<int, int>> screenCoords;
        int minX = INT_MAX, maxX = INT_MIN;
        int minY = INT_MAX, maxY = INT_MIN;
        for (const auto& part : particles) {
            int px = static_cast<int>(viewCenterX + (xCenter - part.x) * scale);
            int py = static_cast<int>(viewCenterY + (yCenter - part.y) * scale);
            screenCoords.emplace_back(px, py);
            minX = std::min(minX, px);
            maxX = std::max(maxX, px);
            minY = std::min(minY, py);
            maxY = std::max(maxY, py);
        }

        // Use a smaller cell size for a higher resolution density grid.
        const int cellSize = 3;
        int gridCols = (maxX - minX) / cellSize + 1;
        int gridRows = (maxY - minY) / cellSize + 1;

        // Create a density grid initialized to zero.
        std::vector<std::vector<int>> density(gridRows, std::vector<int>(gridCols, 0));
        for (const auto& coord : screenCoords) {
            int col = (coord.first - minX) / cellSize;
            int row = (coord.second - minY) / cellSize;
            if (row >= 0 && row < gridRows && col >= 0 && col < gridCols)
                density[row][col]++;
        }

        // Apply a simple 3x3 box blur to spread out intensity.
        std::vector<std::vector<int>> blurredDensity = density; // copy original
        for (int i = 1; i < gridRows - 1; ++i) {
            for (int j = 1; j < gridCols - 1; ++j) {
                int sum = 0;
                for (int di = -1; di <= 1; ++di)
                    for (int dj = -1; dj <= 1; ++dj)
                        sum += density[i + di][j + dj];
                blurredDensity[i][j] = sum / 9; // average value
            }
        }

        // Determine the maximum density from the blurred grid for normalization.
        int maxDensity = 0;
        for (int i = 0; i < gridRows; ++i)
            for (int j = 0; j < gridCols; ++j)
                maxDensity = std::max(maxDensity, blurredDensity[i][j]);

        // Lambda: Map normalized value t (0.0 to 1.0) to a rainbow color.
        auto getRainbowColor = [](double t) -> Color {
            // Define color stops for a full rainbow:
            // 0.0: Purple (128, 0, 128)
            // 0.2: Blue   (0, 0, 255)
            // 0.4: Green  (0, 255, 0)
            // 0.6: Yellow (255, 255, 0)
            // 0.8: Orange (255, 165, 0)
            // 1.0: Red    (255, 0, 0)
            struct ColorStop { double pos; int r, g, b; };
            ColorStop stops[] = {
                {0.0, 0, 0, 0},
                {0.2, 0,   0, 255},
                {0.4, 0, 255,   0},
                {0.6, 255, 255, 0},
                {0.8, 255, 165, 0},
                {1.0, 255, 0,   0}
            };
            if (t < 0.0) t = 0.0;
            if (t > 1.0) t = 1.0;
            ColorStop* left = &stops[0], * right = &stops[1];
            for (int i = 0; i < 5; ++i) {
                if (t >= stops[i].pos && t <= stops[i + 1].pos) {
                    left = &stops[i];
                    right = &stops[i + 1];
                    break;
                }
            }
            double span = right->pos - left->pos;
            double localT = (span > 0.0) ? (t - left->pos) / span : 0.0;
            int r = static_cast<int>(left->r + localT * (right->r - left->r));
            int g = static_cast<int>(left->g + localT * (right->g - left->g));
            int b = static_cast<int>(left->b + localT * (right->b - left->b));
            return Color(255, r, g, b);
            };

        // Create a GDI+ Bitmap for the density grid.
        Bitmap heatmapBmp(gridCols, gridRows, PixelFormat32bppARGB);
        for (int row = 0; row < gridRows; ++row) {
            for (int col = 0; col < gridCols; ++col) {
                double t = (maxDensity > 0) ? (static_cast<double>(blurredDensity[row][col]) / maxDensity) : 0.0;
                Color color = getRainbowColor(t);
                heatmapBmp.SetPixel(col, row, color);
            }
        }

        // Use a Graphics object to scale the heatmap.
        Graphics graphics(hdc);
        graphics.SetSmoothingMode(SmoothingModeHighQuality);
        graphics.SetInterpolationMode(InterpolationModeHighQualityBilinear);

        // Compute the destination rectangle.
        int destWidth = gridCols * cellSize;
        int destHeight = gridRows * cellSize;
        Rect destRect(minX, minY, destWidth, destHeight);

        // Draw the heatmap.
        graphics.DrawImage(&heatmapBmp, destRect);
    }
}


// -----------------------------------------------------------------------------
// Draw Cells: reads "discrete.txt", computes the min/max for normalization,
// and draws each cell as a rectangle with a grayscale color.
// -----------------------------------------------------------------------------
void drawCells(HWND hwnd, HDC hdc, RECT view) {
    std::ifstream file("discrete.txt");
    if (!file)
        return;

    std::vector<std::vector<int>> grid;
    std::string line;
    int datasetCount = 1;
    int globalMin = INT_MAX, globalMax = INT_MIN;
    int totalParticles = 0; // Variable to store the sum

    while (std::getline(file, line)) {
        if (line.empty()) {
            if (++datasetCount > dataset)
                break;
            continue;
        }
        if (datasetCount == dataset) {
            std::istringstream iss(line);
            std::vector<int> row;
            int value;
            while (iss >> value) {
                row.push_back(value);
                if (value > 0) {
                    globalMin = std::min(globalMin, value);
                    globalMax = std::max(globalMax, value);
                }
                totalParticles += value; // Accumulate the sum
            }
            grid.push_back(row);
        }
    }
    file.close();
    if (grid.empty())
        return;

    RECT clientRect;
    GetClientRect(hwnd, &clientRect);

    int numRows = grid.size();
    int numCols = grid[0].size();
    int viewCenterX = (view.right + view.left) / 2;
    int viewCenterY = (view.bottom + view.top) / 2;

    int cellSize = static_cast<int>((clientRect.bottom - clientRect.top) / static_cast<double>(numCols) * 0.65);
    int range = std::max(1, globalMax - globalMin);

    int lineLength = cellSize;

    int gridLeft = viewCenterX - (numRows * cellSize) / 2;
    int gridTop = viewCenterY - (numCols * cellSize) / 2;
    int gridRight = gridLeft + numRows * cellSize;
    int gridBottom = gridTop + numCols * cellSize;

    HPEN pen = CreatePen(PS_SOLID, 1, RGB(255, 255, 255));
    HPEN oldPen = (HPEN)SelectObject(hdc, pen);

    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            int newI = numCols - 1 - j;
            int newJ = numRows - 1 - i;
            int left = viewCenterX - (numRows * cellSize) / 2 + newJ * cellSize;
            int top = viewCenterY - (numCols * cellSize) / 2 + newI * cellSize;
            int gray = (grid[i][j] == 0) ? 0 : static_cast<int>(255 * (grid[i][j] - globalMin) / static_cast<double>(range));
            COLORREF color = (grid[i][j] == 0) ? RGB(0, 0, 0) : RGB(gray, gray, gray);
            HBRUSH brush = CreateSolidBrush(color);
            RECT cellRect = { left, top, left + cellSize, top + cellSize };
            FillRect(hdc, &cellRect, brush);
            DeleteObject(brush);

            // Draw Grid if options is selected
            if (showGridPoints)
            {
                // Draw small corner lines
                MoveToEx(hdc, left, top, nullptr); // Top-left
                LineTo(hdc, left + lineLength / 6, top);
                MoveToEx(hdc, left, top, nullptr);
                LineTo(hdc, left, top + lineLength / 6);

                MoveToEx(hdc, left + cellSize, top, nullptr); // Top-right
                LineTo(hdc, left + cellSize - lineLength / 6, top);
                MoveToEx(hdc, left + cellSize, top, nullptr);
                LineTo(hdc, left + cellSize, top + lineLength / 6);

                MoveToEx(hdc, left, top + cellSize, nullptr); // Bottom-left
                LineTo(hdc, left + lineLength / 6, top + cellSize);
                MoveToEx(hdc, left, top + cellSize, nullptr);
                LineTo(hdc, left, top + cellSize - lineLength / 6);

                MoveToEx(hdc, left + cellSize, top + cellSize, nullptr); // Bottom-right
                LineTo(hdc, left + cellSize - lineLength / 6, top + cellSize);
                MoveToEx(hdc, left + cellSize, top + cellSize, nullptr);
                LineTo(hdc, left + cellSize, top + cellSize - lineLength / 6);
            }

            // Draws grid corners
            int cornerSize = 7;
            MoveToEx(hdc, gridLeft, gridBottom, nullptr); // Bottom-Left
            LineTo(hdc, gridLeft + 7, gridBottom);
            MoveToEx(hdc, gridLeft, gridBottom, nullptr);
            LineTo(hdc, gridLeft, gridBottom - 7);
            MoveToEx(hdc, gridRight, gridBottom, nullptr); // Bottom-Right
            LineTo(hdc, gridRight - 7, gridBottom);
            MoveToEx(hdc, gridRight, gridBottom, nullptr);
            LineTo(hdc, gridRight, gridBottom - 7);
            MoveToEx(hdc, gridRight, gridTop, nullptr); // Top-Right
            LineTo(hdc, gridRight - 7, gridTop);
            MoveToEx(hdc, gridRight, gridTop, nullptr);
            LineTo(hdc, gridRight, gridTop + 7);
            MoveToEx(hdc, gridLeft, gridTop, nullptr); // Top-Left
            LineTo(hdc, gridLeft + 7, gridTop);
            MoveToEx(hdc, gridLeft, gridTop, nullptr);
            LineTo(hdc, gridLeft, gridTop + 7);
        }
    }

    SelectObject(hdc, oldPen);
    DeleteObject(pen);

    std::wstring datasetText = L"DtsN: " + std::to_wstring(dataset);
    RECT textRect = { view.right - 210, view.bottom - 30, view.right - 10, view.bottom - 10 };
    SetTextColor(hdc, RGB(255, 255, 255));
    SetBkMode(hdc, TRANSPARENT);
    DrawTextW(hdc, datasetText.c_str(), -1, &textRect, DT_RIGHT | DT_SINGLELINE);

    std::wostringstream woss;
    woss << L"NoP: " << std::scientific << std::setprecision(2) << static_cast<double>(totalParticles);
    std::wstring particleText = woss.str();
    RECT particleTextRect = { view.right - 210, view.bottom - 50, view.right - 10, view.bottom - 30 };
    DrawTextW(hdc, particleText.c_str(), -1, &particleTextRect, DT_RIGHT | DT_SINGLELINE);
}


// -----------------------------------------------------------------------------
// Draw UI Controls: buttons and radio inputs using GDI+ for anti-aliasing.
// -----------------------------------------------------------------------------
void drawRadioButton(HDC hdc, int centerX, int centerY, int radius, bool filled) {
    Graphics graphics(hdc);
    graphics.SetSmoothingMode(SmoothingModeAntiAlias);
    Pen pen(Color(255, 0, 0, 0), 2);
    SolidBrush brush(filled ? Color(255, 0, 0, 0) : Color(255, 200, 200, 200));
    graphics.FillEllipse(&brush, centerX - radius / 2, centerY - radius / 2, radius, radius);
    graphics.DrawEllipse(&pen, centerX - radius / 2, centerY - radius / 2, radius, radius);
    graphics.Flush();  // Ensure the drawing is flushed to the DC
}

void drawButtonsAndInputs(HDC hdc, const RECT& view3) {
    int buffer = BUTTON_BUFFER;
    int bWidth = BUTTON_WIDTH;
    int buttonHeight = (view3.bottom - view3.top) - 2 * buffer;
    int centerY = (view3.top + view3.bottom) / 2;

    // Draw Particles button
    {
        int centerX = buffer + bWidth / 2;
        RECT btnRect = { centerX - bWidth / 2, centerY - buttonHeight / 2, centerX + bWidth / 2, centerY + buttonHeight / 2 };
        HBRUSH brush = CreateSolidBrush(RGB(200, 200, 200));
        FillRect(hdc, &btnRect, brush);
        SetTextColor(hdc, RGB(0, 0, 0));
        SetBkMode(hdc, TRANSPARENT);
        DrawText(hdc, "Particles (C)", -1, &btnRect, DT_CENTER | DT_VCENTER | DT_SINGLELINE);
        DeleteObject(brush);
    }
    // Draw Cells button
    {
        int centerX = buffer + bWidth + buffer + bWidth / 2;
        RECT btnRect = { centerX - bWidth / 2, centerY - buttonHeight / 2, centerX + bWidth / 2, centerY + buttonHeight / 2 };
        HBRUSH brush = CreateSolidBrush(RGB(200, 200, 200));
        FillRect(hdc, &btnRect, brush);
        SetTextColor(hdc, RGB(0, 0, 0));
        SetBkMode(hdc, TRANSPARENT);
        DrawText(hdc, "Cells (Step)", -1, &btnRect, DT_CENTER | DT_VCENTER | DT_SINGLELINE);
        DeleteObject(brush);
    }
    // If in cells mode, draw radio buttons and scrollbar
    if (drawCellsFlag) {
        drawScrollBar(hdc, computedScrollBarX, computedScrollBarY);
    }
    if (drawParticlesFlag) {
        int radioX = 2 * bWidth + 3 * buffer;
        int radioY1 = centerY - BUTTON_BUFFER;
        int radioY2 = centerY + BUTTON_BUFFER;
        drawRadioButton(hdc, radioX, radioY1, 10, optionsCells == 1);
        drawRadioButton(hdc, radioX, radioY2, 10, optionsCells == 2);
    }
}

// -----------------------------------------------------------------------------
// Window Procedure: handles resizing, painting, mouse input, etc.
// -----------------------------------------------------------------------------
LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam) {
    static RECT view2, view3, view4;
    static int viewCenterX = 0, viewCenterY = 0;
    int width = LOWORD(lParam);
    int height = HIWORD(lParam);

    switch (uMsg)
    {
    case WM_SIZE:
    {
        const int border1 = 5, border2 = 8, border3 = 32;
        view2 = { 0, height / border2, width, height };
        view3 = { 0, height / border3, width, height / border2 };
        view4 = { 0, 0, width, height / border3 };
        viewCenterX = (view2.left + view2.right) / 2;
        viewCenterY = (view2.top + view2.bottom) / 2;
        computedScrollBarY = view2.bottom - SCROLL_BAR_HEIGHT - MARGIN_BOTTOM;
        computedScrollBarX = view2.left + ((width - SCROLL_BAR_WIDTH) / 2);
        InvalidateRect(hwnd, NULL, TRUE);
        return 0;
    }
    case WM_LBUTTONDOWN:
    {
        int mouseX = GET_X_LPARAM(lParam);
        int mouseY = GET_Y_LPARAM(lParam);
        int buffer = BUTTON_BUFFER;
        int bWidth = BUTTON_WIDTH;
        int buttonTop = view3.top + buffer;
        int buttonBottom = buttonTop + ((view3.bottom - view3.top) - 2 * buffer);
        // Detect click in Particles button area.
        if (drawCellsFlag && mouseX >= buffer && mouseX <= buffer + bWidth &&
            mouseY >= buttonTop && mouseY <= buttonBottom) {
            drawParticlesFlag = true;
            drawCellsFlag = false;
            InvalidateRect(hwnd, NULL, TRUE);
        }
        // Detect click in Cells button area.
        if (drawParticlesFlag && mouseX >= (buffer + bWidth + buffer) && mouseX <= (buffer + bWidth + buffer + bWidth) &&
            mouseY >= buttonTop && mouseY <= buttonBottom) {
            drawCellsFlag = true;
            drawParticlesFlag = false;
            InvalidateRect(hwnd, NULL, TRUE);
        }
        if (drawParticlesFlag)
        {
            int bWidth = BUTTON_WIDTH;
            int buffer = BUTTON_BUFFER;
            // Calculate radioX the same way as in drawButtonsAndInputs.
            int radioX = 2 * bWidth + 3 * buffer;
            // Calculate centerY from view3.
            int centerY = (view3.top + view3.bottom) / 2;
            // Define radio button centers (offset vertically by BUTTON_BUFFER).
            int radioY1 = centerY - BUTTON_BUFFER;
            int radioY2 = centerY + BUTTON_BUFFER;

            // Check if the click is within a 10-pixel radius of the first radio button.
            if (std::abs(mouseX - radioX) <= 10 && std::abs(mouseY - radioY1) <= 10)
            {
                radio = 1;
            }
            // Otherwise, check the second radio button.
            else if (std::abs(mouseX - radioX) <= 10 && std::abs(mouseY - radioY2) <= 10)
            {
                radio = 2;
            }
            InvalidateRect(hwnd, NULL, TRUE);
        }
        // Check for scrollbar thumb dragging.
        int thumbPos = computedScrollBarX;
        if (totalDatasets > 1)
            thumbPos += (scrollValue * (SCROLL_BAR_WIDTH - THUMB_WIDTH)) / (totalDatasets - 1);
        if (drawCellsFlag && mouseX >= thumbPos && mouseX <= thumbPos + THUMB_WIDTH &&
            mouseY >= computedScrollBarY && mouseY <= computedScrollBarY + SCROLL_BAR_HEIGHT) {
            isDragging = true;
            dragStartX = mouseX;
            thumbStartPos = thumbPos;
            SetCapture(hwnd);
        }
        return 0;
    }
    case WM_MOUSEMOVE:
    {
        if (isDragging) {
            int mouseX = GET_X_LPARAM(lParam);
            int deltaX = mouseX - dragStartX;
            int newThumbPos = thumbStartPos + deltaX;
            newThumbPos = std::clamp(newThumbPos, computedScrollBarX, computedScrollBarX + SCROLL_BAR_WIDTH - THUMB_WIDTH);
            if (totalDatasets > 1)
                scrollValue = ((newThumbPos - computedScrollBarX) * (totalDatasets - 1)) / (SCROLL_BAR_WIDTH - THUMB_WIDTH);
            InvalidateRect(hwnd, NULL, TRUE);
        }
        return 0;
    }
    case WM_LBUTTONUP:
    {
        if (isDragging) {
            isDragging = false;
            ReleaseCapture();
        }
        return 0;
    }
    case WM_PAINT:
    {
        PAINTSTRUCT ps;
        HDC hdc = BeginPaint(hwnd, &ps);
        RECT clientRect;
        GetClientRect(hwnd, &clientRect);
        int cWidth = clientRect.right - clientRect.left;
        int cHeight = clientRect.bottom - clientRect.top;

        // Create a memory device context for double-buffering.
        HDC memDC = CreateCompatibleDC(hdc);
        HBITMAP memBitmap = CreateCompatibleBitmap(hdc, cWidth, cHeight);
        HGDIOBJ oldBitmap = SelectObject(memDC, memBitmap);

        // Clear background.
        FillRect(memDC, &clientRect, (HBRUSH)GetStockObject(BLACK_BRUSH));

        // Draw view sections.
        HBRUSH brushGray = CreateSolidBrush(RGB(0, 0, 0));
        HBRUSH brushSet = CreateSolidBrush(RGB(80, 80, 80));
        HBRUSH brushBar = CreateSolidBrush(RGB(40, 40, 40));
        FillRect(memDC, &view2, brushGray);
        FillRect(memDC, &view3, brushSet);
        FillRect(memDC, &view4, brushBar);
        DeleteObject(brushGray);
        DeleteObject(brushSet);
        DeleteObject(brushBar);

        // Draw buttons and additional UI elements.
        drawButtonsAndInputs(memDC, view3);

        // Draw particles or cells based on current mode.
        if (drawParticlesFlag)
            drawParticles(hwnd, memDC, viewCenterX, viewCenterY);
        if (drawCellsFlag)
            drawCells(hwnd, memDC, view2);

        // Blit the memory DC to screen.
        BitBlt(hdc, 0, 0, cWidth, cHeight, memDC, 0, 0, SRCCOPY);
        SelectObject(memDC, oldBitmap);
        DeleteObject(memBitmap);
        DeleteDC(memDC);
        EndPaint(hwnd, &ps);
        return 0;
    }
    case WM_DESTROY:
        ShutdownGDIPlus();
        PostQuitMessage(0);
        return 0;
    default:
        return DefWindowProc(hwnd, uMsg, wParam, lParam);
    }
}

// -----------------------------------------------------------------------------
// WinMain: Application entry point.
// -----------------------------------------------------------------------------
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE, LPSTR, int nCmdShow)
{
    const char CLASS_NAME[] = "ImprovedWindowClass";
    WNDCLASSEX wc = { 0 };
    wc.cbSize = sizeof(wc);
    wc.lpfnWndProc = WindowProc;
    wc.hInstance = hInstance;
    wc.lpszClassName = CLASS_NAME;
    if (!RegisterClassEx(&wc))
        return 0;

    HWND hwnd = CreateWindowEx(
        0,
        CLASS_NAME,
        "DCFHT (PIC SIM)",
        WS_OVERLAPPEDWINDOW,
        CW_USEDEFAULT, CW_USEDEFAULT, 800, 600,
        nullptr, nullptr, hInstance, nullptr
    );
    if (!hwnd)
        return 0;

    HICON hIcon = (HICON)LoadImage(NULL, "icon.ico", IMAGE_ICON, 32, 32, LR_LOADFROMFILE);
    if (hIcon) {
        SendMessage(hwnd, WM_SETICON, ICON_BIG, (LPARAM)hIcon);
        SendMessage(hwnd, WM_SETICON, ICON_SMALL, (LPARAM)hIcon);
    }

    ShowWindow(hwnd, nCmdShow);
    UpdateWindow(hwnd);
    InitGDIPlus();

    // Update the totalDatasets value using the countDatasets function.
    totalDatasets = countDatasets("discrete.txt");
    if (totalDatasets < 1)
        totalDatasets = 1; // Ensure there's at least one dataset

    MSG msg;
    while (GetMessage(&msg, nullptr, 0, 0) > 0) {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
    return 0;
}
