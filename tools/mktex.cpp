// mktex — CC0 のネオンテトラ写真から「展開したボディテクスチャ」を作る。
//
// 写真の魚は体高が場所によって違うので、そのまま貼ると伸び縮みする。
// そこで体の上縁・下縁を制御点で与え、各列を上縁〜下縁で正規化して長方形に開く(unwrap)。
// こうしておけば、シミュレータ側は s(鼻先→尾柄) と v(背-1〜腹+1) で素直に引ける。
//
//   clang++ -O2 -std=c++17 -I. tools/mktex.cpp -o /tmp/mktex
//   /tmp/mktex assets/neon_tetra_cc0.jpg src/tetra_tex.h /tmp/unwrap.png
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>

static const int TW = 224, TH = 72;    // 出力テクスチャ

// 元画像(2372x1476)での体の輪郭。魚は右を向いている。
// x が大きいほど鼻先。fin は含めず「胴体」だけ。
struct P { double x, y; };
static const P TOP[] = {
    {2077, 524}, {1975, 417}, {1821, 375}, {1583, 339}, {1345, 321},
    {1142, 345}, { 964, 393}, { 797, 470}, { 666, 547}, { 559, 619},
    { 506, 666}
};
static const P BOT[] = {
    {2077, 541}, {1975, 619}, {1821, 678}, {1583, 732}, {1345, 768},
    {1142, 785}, { 964, 800}, { 797, 809}, { 666, 812}, { 559, 807},
    { 506, 800}
};
static const int NP = (int)(sizeof(TOP) / sizeof(TOP[0]));
// 輪郭ぴったりだと背景を拾うので、少し内側に詰めてから展開する
static const double INSET = 0.06;

// 制御点列を x で線形補間
static double edge_y(const P* pts, double x) {
    if (x >= pts[0].x) return pts[0].y;
    for (int i = 0; i < NP - 1; ++i) {
        if (x <= pts[i].x && x >= pts[i + 1].x) {
            double t = (pts[i].x - x) / (pts[i].x - pts[i + 1].x);
            t = t * t * (3 - 2 * t);
            return pts[i].y * (1 - t) + pts[i + 1].y * t;
        }
    }
    return pts[NP - 1].y;
}

int main(int argc, char** argv) {
    const char* in  = argc > 1 ? argv[1] : "assets/neon_tetra_cc0.jpg";
    const char* out = argc > 2 ? argv[2] : "src/tetra_tex.h";
    const char* dbg = argc > 3 ? argv[3] : nullptr;

    int w, h, n;
    unsigned char* img = stbi_load(in, &w, &h, &n, 3);
    if (!img) { fprintf(stderr, "load failed: %s\n", in); return 1; }
    fprintf(stderr, "loaded %dx%d\n", w, h);

    std::vector<unsigned char> tex((size_t)TW * TH * 4, 0);
    double x0 = TOP[NP - 1].x, x1 = TOP[0].x;          // 尾柄 → 鼻先

    for (int i = 0; i < TW; ++i) {
        double s = (i + 0.5) / TW;                      // 0=鼻先, 1=尾柄
        double sx = x1 + (x0 - x1) * s;
        double yt0 = edge_y(TOP, sx), yb0 = edge_y(BOT, sx);
        double mid = (yt0 + yb0) * 0.5, hh2 = (yb0 - yt0) * 0.5 * (1.0 - INSET);
        double yt = mid - hh2, yb = mid + hh2;
        for (int j = 0; j < TH; ++j) {
            double v = (j + 0.5) / TH;                  // 0=背, 1=腹
            // 上下に少しはみ出して縁をなだらかに取る
            double sy = yt + (yb - yt) * v;
            int px = (int)(sx + 0.5), py = (int)(sy + 0.5);
            if (px < 0) px = 0; if (px >= w) px = w - 1;
            if (py < 0) py = 0; if (py >= h) py = h - 1;
            const unsigned char* p = &img[((size_t)py * w + px) * 3];
            size_t o = ((size_t)j * TW + i) * 4;
            tex[o + 0] = p[0]; tex[o + 1] = p[1]; tex[o + 2] = p[2]; tex[o + 3] = 255;
        }
    }
    stbi_image_free(img);

    if (dbg) {
        // 展開結果を4倍に拡大して確認用に出す
        int dw = TW * 4, dh = TH * 4;
        std::vector<unsigned char> big((size_t)dw * dh * 4);
        for (int j = 0; j < dh; ++j) for (int i = 0; i < dw; ++i) {
            size_t o = ((size_t)j * dw + i) * 4, q = ((size_t)(j / 4) * TW + (i / 4)) * 4;
            for (int c = 0; c < 4; ++c) big[o + c] = tex[q + c];
        }
        stbi_write_png(dbg, dw, dh, 4, big.data(), dw * 4);
        // 元写真に輪郭を重ねた確認用画像
        int w2, h2, n2;
        unsigned char* im2 = stbi_load(in, &w2, &h2, &n2, 3);
        if (im2) {
            std::vector<unsigned char> ov((size_t)w2 * h2 * 3);
            for (size_t k = 0; k < ov.size(); ++k) ov[k] = im2[k];
            for (double sx2 = x0; sx2 <= x1; sx2 += 0.5) {
                double yt = edge_y(TOP, sx2), yb = edge_y(BOT, sx2);
                for (int d = -2; d <= 2; ++d) {
                    int xi = (int)sx2;
                    int a = (int)yt + d, b = (int)yb + d;
                    if (xi >= 0 && xi < w2 && a >= 0 && a < h2) {
                        size_t o = ((size_t)a * w2 + xi) * 3; ov[o] = 0; ov[o+1] = 255; ov[o+2] = 0; }
                    if (xi >= 0 && xi < w2 && b >= 0 && b < h2) {
                        size_t o = ((size_t)b * w2 + xi) * 3; ov[o] = 255; ov[o+1] = 0; ov[o+2] = 255; }
                }
            }
            stbi_write_png("/tmp/tetra_outline.png", w2, h2, 3, ov.data(), w2 * 3);
            stbi_image_free(im2);
        }
    }

    FILE* f = fopen(out, "w");
    if (!f) { fprintf(stderr, "cannot write %s\n", out); return 1; }
    fprintf(f, "// 自動生成 (tools/mktex.cpp)。CC0 のネオンテトラ写真を胴体で展開したテクスチャ。\n");
    fprintf(f, "// 出典: https://www.goodfreephotos.com/animals/fish/neon-tetra.jpg.php (CC0 / Public Domain)\n");
    fprintf(f, "static const int TETRA_TW = %d, TETRA_TH = %d;\n", TW, TH);
    fprintf(f, "static const unsigned char TETRA_TEX[] = {\n");
    for (size_t k = 0; k < tex.size(); ++k) {
        fprintf(f, "%u,", tex[k]);
        if ((k % 32) == 31) fprintf(f, "\n");
    }
    fprintf(f, "};\n");
    fclose(f);
    fprintf(stderr, "wrote %s (%dx%d)\n", out, TW, TH);
    return 0;
}
