// mktex — 金魚の写真から「展開したボディテクスチャ」を作る。
//
// 写真の魚は場所によって体高が違うので、そのまま貼ると伸び縮みする。
// そこで**胴体の上縁・下縁を自動で拾い**、各列を上縁〜下縁で正規化して長方形に開く(unwrap)。
// こうしておけば、シミュレータ側は s(鼻先0→尾柄1) と v(背-1〜腹+1) で素直に引ける。
//
// 輪郭を手で測って制御点で与える作りだと、魚を変えるたびに測り直しになる。
// 出目金やピンポンパールも同じ道具で通したいので、自動にしてある。分け方は2段:
//
//   1) 魚と背景 … α があればそれを使う。無ければ「四隅の色に近いもの」を背景とみなす
//   2) 胴体とひれ … **ひれは胴体より明るい**(金魚の胴は濃い朱、ひれは淡い黄)。
//      明るさのヒストグラムを大津の判別分析で2つに割り、暗いほうを胴体とする。
//      しきい値を決め打ちにすると写真ごとに破綻するので、写真から決めさせる
//
//   clang++ -O2 -std=c++17 -I. tools/mktex.cpp -o /tmp/mktex
//   /tmp/mktex kingyo.png src/kingyo_tex.h /tmp/unwrap.png /tmp/outline.png KINGYO
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>


int main(int argc, char** argv) {
    const char* in   = argc > 1 ? argv[1] : "kingyo.png";
    const char* out  = argc > 2 ? argv[2] : "src/kingyo_tex.h";
    const char* dbg  = argc > 3 ? argv[3] : nullptr;   // 作ったスプライトの確認用
    const char* dbg2 = argc > 4 ? argv[4] : nullptr;   // 元写真＋拾った輪郭
    const char* sym  = argc > 5 ? argv[5] : "KINGYO";

    int w, h, n;
    unsigned char* img = stbi_load(in, &w, &h, &n, 4);   // α も読む
    if (!img) { fprintf(stderr, "load failed: %s\n", in); return 1; }
    fprintf(stderr, "loaded %dx%d (%dch)\n", w, h, n);
    auto A = [&](int x, int y, int c) { return (int)img[((size_t)y * w + x) * 4 + c]; };

    // --- 1) 魚と背景
    std::vector<unsigned char> fishm((size_t)w * h, 0);
    bool hasAlpha = false;
    for (size_t i = 0; i < (size_t)w * h && !hasAlpha; ++i) if (img[i * 4 + 3] < 250) hasAlpha = true;
    int bgr = 0, bgg = 0, bgb = 0;
    if (!hasAlpha) {    // 四隅の平均を背景色とみなす
        int cx[4] = {0, w - 1, 0, w - 1}, cy[4] = {0, 0, h - 1, h - 1};
        for (int k = 0; k < 4; ++k) { bgr += A(cx[k],cy[k],0); bgg += A(cx[k],cy[k],1); bgb += A(cx[k],cy[k],2); }
        bgr /= 4; bgg /= 4; bgb /= 4;
    }
    long nfish = 0;
    for (int y = 0; y < h; ++y) for (int x = 0; x < w; ++x) {
        bool on;
        if (hasAlpha) on = A(x, y, 3) > 128;
        else {
            int dr = A(x,y,0) - bgr, dg = A(x,y,1) - bgg, db = A(x,y,2) - bgb;
            on = (dr*dr + dg*dg + db*db) > 42 * 42;
        }
        fishm[(size_t)y * w + x] = on ? 1 : 0;
        if (on) nfish++;
    }
    fprintf(stderr, "魚の画素 %ld (%s)\n", nfish, hasAlpha ? "α から" : "背景色との差から");
    if (nfish < 500) { fprintf(stderr, "魚が見つからない\n"); return 1; }

    // --- 2) 胴体とひれ。明るさのヒストグラムを大津の判別分析で割る
    int hist[256] = {0};
    for (int y = 0; y < h; ++y) for (int x = 0; x < w; ++x)
        if (fishm[(size_t)y * w + x]) hist[(A(x,y,0) * 30 + A(x,y,1) * 59 + A(x,y,2) * 11) / 100]++;
    long tot = nfish, sum = 0;
    for (int i = 0; i < 256; ++i) sum += (long)i * hist[i];
    long wB = 0, sumB = 0; double best = -1; int thr = 128;
    for (int i = 0; i < 256; ++i) {
        wB += hist[i]; if (!wB) continue;
        long wF = tot - wB; if (!wF) break;
        sumB += (long)i * hist[i];
        double mB = (double)sumB / wB, mF = (double)(sum - sumB) / wF;
        double v = (double)wB * wF * (mB - mF) * (mB - mF);
        if (v > best) { best = v; thr = i; }
    }
    fprintf(stderr, "胴体とひれの境目(明るさ) = %d\n", thr);
    std::vector<unsigned char> bodym((size_t)w * h, 0);
    for (int y = 0; y < h; ++y) for (int x = 0; x < w; ++x)
        if (fishm[(size_t)y * w + x] &&
            (A(x,y,0) * 30 + A(x,y,1) * 59 + A(x,y,2) * 11) / 100 <= thr)
            bodym[(size_t)y * w + x] = 1;

    // --- 列ごとに「いちばん長く続いている胴体の縦の並び」を取る。
    //     ひれは胴体から離れた別の並びになるので、これで落ちる
    std::vector<int> topY(w, -1), botY(w, -1), runL(w, 0);
    for (int x = 0; x < w; ++x) {
        int bs = -1, bl = 0, cs = -1, cl = 0;
        for (int y = 0; y <= h; ++y) {
            int on = (y < h) ? bodym[(size_t)y * w + x] : 0;
            if (on) { if (cs < 0) cs = y; cl++; }
            else { if (cl > bl) { bl = cl; bs = cs; } cs = -1; cl = 0; }
        }
        if (bl > 0) { topY[x] = bs; botY[x] = bs + bl - 1; runL[x] = bl; }
    }
    int maxRun = 0; for (int x = 0; x < w; ++x) maxRun = std::max(maxRun, runL[x]);
    // 胴体の x 範囲 = 十分な高さがある列が続いているところ(いちばん長い区間)。
    // しきい値を上げすぎると**鼻先や尾柄の細いところが落ちる**ので低めに取る
    int bestA = -1, bestB = -1, curA = -1;
    for (int x = 0; x <= w; ++x) {
        bool ok = (x < w) && runL[x] > maxRun * 0.10;
        if (ok) { if (curA < 0) curA = x; }
        else { if (curA >= 0 && (bestA < 0 || x - curA > bestB - bestA)) { bestA = curA; bestB = x; } curA = -1; }
    }
    if (bestA < 0) { fprintf(stderr, "胴体が見つからない\n"); return 1; }
    fprintf(stderr, "胴体の x 範囲 %d〜%d (最大の体高 %d)\n", bestA, bestB - 1, maxRun);

    // 縁をならす。**ひれの付け根で輪郭が凹む**のをここで埋める。
    // 背びれ・尻びれの根元は明るいので胴体から外れ、その列だけ縁が体の内側に食い込む。
    // 広い窓でならした線(trend)より内側に入っている列は、trend で置き換える。
    // 単純に平均するだけだと凹みごと均されて、体が痩せてしまう
    std::vector<double> tS(w, 0), bS(w, 0);
    auto smooth = [&](const std::vector<double>& src, std::vector<double>& dst, int win) {
        for (int x = bestA; x < bestB; ++x) {
            double sum = 0; int c = 0;
            for (int k = -win; k <= win; ++k) {
                int xx = std::min(bestB - 1, std::max(bestA, x + k));
                sum += src[xx]; c++;
            }
            dst[x] = sum / c;
        }
    };
    std::vector<double> tRaw(w, 0), bRaw(w, 0), tTr(w, 0), bTr(w, 0);
    for (int x = bestA; x < bestB; ++x) { tRaw[x] = topY[x]; bRaw[x] = botY[x]; }
    const int WID = std::max(4, (bestB - bestA) / 6);      // 広い窓 = ひれの基底より広く
    for (int it = 0; it < 3; ++it) {
        smooth(tRaw, tTr, WID); smooth(bRaw, bTr, WID);
        for (int x = bestA; x < bestB; ++x) {
            if (tRaw[x] > tTr[x]) tRaw[x] = tTr[x];        // 上へ食い込んでいたら戻す
            if (bRaw[x] < bTr[x]) bRaw[x] = bTr[x];        // 下も同じ
        }
    }
    smooth(tRaw, tS, std::max(2, (bestB - bestA) / 40));   // 最後に軽くならす
    smooth(bRaw, bS, std::max(2, (bestB - bestA) / 40));
    // どちらが頭か。**端の平均で比べてはいけない**: 鼻先も尾柄も細いので差が出ない。
    // 金魚は頭側だけ内側へ入るとすぐ深くなる(肩)。少し内側(2割)の体高で比べる
    int q = std::max(1, (bestB - bestA) / 5);
    double hL = bS[bestA + q] - tS[bestA + q];      // 左端から2割の体高
    double hR = bS[bestB - 1 - q] - tS[bestB - 1 - q];
    bool headRight = hR > hL;
    fprintf(stderr, "体高: 左から2割 %.0f / 右から2割 %.0f\n", hL, hR);
    fprintf(stderr, "頭は %s向き\n", headRight ? "右" : "左");

    if (dbg2) {   // 元写真に拾った輪郭を重ねる
        std::vector<unsigned char> ov((size_t)w * h * 3, 0);
        for (int y = 0; y < h; ++y) for (int x = 0; x < w; ++x) {
            size_t o = ((size_t)y * w + x) * 3;
            for (int c = 0; c < 3; ++c) ov[o + c] = (unsigned char)A(x, y, c);
            if (!fishm[(size_t)y * w + x]) { ov[o] = 30; ov[o+1] = 30; ov[o+2] = 40; }
        }
        for (int x = bestA; x < bestB; ++x) {
            for (int k = -1; k <= 1; ++k) {
                int yt = (int)(tS[x] + 0.5) + k, yb = (int)(bS[x] + 0.5) + k;
                if (yt >= 0 && yt < h) { size_t o = ((size_t)yt * w + x) * 3; ov[o]=0; ov[o+1]=255; ov[o+2]=90; }
                if (yb >= 0 && yb < h) { size_t o = ((size_t)yb * w + x) * 3; ov[o]=0; ov[o+1]=200; ov[o+2]=255; }
            }
        }
        stbi_write_png(dbg2, w, h, 3, ov.data(), w * 3);
    }

    const double bodyLenGlobal = (double)(bestB - bestA);
    // ---- 体の輪郭も書き出す。**手で決めた体型に写真のひれを付けても、付け根が合わない。**
    // 中心線は「鼻先の中点〜尾柄の中点」を結んだ直線。そこからの上下の張り出しを体長で割る。
    // 魚は背と腹で膨らみ方が違うので、上下は別々に持つ(左右対称の体型にすると別の魚になる)
    const int NSEG = 32;
    std::vector<double> profT(NSEG + 1), profB(NSEG + 1);
    {
        int xn = headRight ? (bestB - 1) : bestA;
        int xp = headRight ? bestA : (bestB - 1);
        double yn = (tS[xn] + bS[xn]) * 0.5, yp = (tS[xp] + bS[xp]) * 0.5;
        for (int i2 = 0; i2 <= NSEG; ++i2) {
            double u = (double)i2 / NSEG;
            double xf = xn + (xp - xn) * u;
            int x = std::min(bestB - 1, std::max(bestA, (int)(xf + 0.5)));
            double axisY = yn + (yp - yn) * u;
            // 展開テクスチャは「上縁と下縁の中点」を v=0 にしているので、
            // 貼る側もそこを基準にする。中心線(鼻先〜尾柄の直線)からのずれと、半分の高さ
            profT[i2] = ((tS[x] + bS[x]) * 0.5 - axisY) / bodyLenGlobal;   // 中点のずれ
            profB[i2] = (bS[x] - tS[x]) * 0.5 / bodyLenGlobal;             // 半分の高さ
        }
        double mt = 0, mb = 0;
        for (int i2 = 0; i2 <= NSEG; ++i2) { mt = std::max(mt, profT[i2]); mb = std::max(mb, profB[i2]); }
        fprintf(stderr, "体の輪郭: 張り出しの最大 背 %.3f / 腹 %.3f (体長比)\n", mt, mb);
    }

    // ---- 魚まるごと1枚のスプライト。
    // 体を短冊に展開し、ひれを別に貼る作りは、写真の「魚らしさ」を分解してしまう。
    // 頭の丸み、ひれの付き方、目の位置がそれぞれ別々に近似され、どこかが必ず嘘になる。
    // 写真を1枚のまま格子で曲げれば、頭もひれも目も写真のまま動く。
    const int SW = 448, SH = 224;
    std::vector<unsigned char> full((size_t)SW * SH * 4, 0);
    double sprNose = 0, sprPed = 1, sprAxis = 0.5, sprScale = 1, sprAspect = 0.5;
    {
        int ax2 = w, bx2 = -1, ay2 = h, by2 = -1;
        for (int y = 0; y < h; ++y) for (int x = 0; x < w; ++x)
            if (fishm[(size_t)y * w + x]) { ax2 = std::min(ax2,x); bx2 = std::max(bx2,x);
                                            ay2 = std::min(ay2,y); by2 = std::max(by2,y); }
        double fx0 = ax2, fx1 = bx2, fy0 = ay2, fy1 = by2;
        for (int j = 0; j < SH; ++j) for (int i = 0; i < SW; ++i) {
            // 頭が左に来るように詰める(貼る側は「鼻先→尾」の向きだけ考えればよくなる)
            double u = (double)i / (SW - 1), v = (double)j / (SH - 1);
            double sx = headRight ? (fx1 - u * (fx1 - fx0)) : (fx0 + u * (fx1 - fx0));
            double sy = fy0 + v * (fy1 - fy0);
            int x = std::min(w-1, std::max(0, (int)(sx + 0.5)));
            int y = std::min(h-1, std::max(0, (int)(sy + 0.5)));
            size_t o = (size_t)(j * SW + i) * 4;
            full[o+0] = (unsigned char)A(x,y,0);
            full[o+1] = (unsigned char)A(x,y,1);
            full[o+2] = (unsigned char)A(x,y,2);
            // ひれの先は薄い。白い紙の上に透けているとみて α を見積もり、
            // **背景の白を割り戻す**。割り戻さないと、ひれが「白っぽい膜」になり、
            // 縁に白いふちが残る (obs = a*C + (1-a)*255 を C について解く)
            double a2;
            if (!fishm[(size_t)y * w + x]) a2 = 0.0;
            else if (hasAlpha) a2 = A(x,y,3) / 255.0;
            else {
                double mn = std::min(A(x,y,0), std::min(A(x,y,1), A(x,y,2)));
                double raw = 1.0 - mn / 255.0;
                if (raw < 0.05) a2 = 0.0;                     // ほぼ白 = 背景の写り込み
                else {
                    a2 = std::min(1.0, raw * 1.5 + 0.15);
                    double inv = 1.0 / a2, wht = (1.0 - a2) * 255.0;
                    full[o+0] = (unsigned char)std::min(255.0, std::max(0.0, (A(x,y,0) - wht) * inv));
                    full[o+1] = (unsigned char)std::min(255.0, std::max(0.0, (A(x,y,1) - wht) * inv));
                    full[o+2] = (unsigned char)std::min(255.0, std::max(0.0, (A(x,y,2) - wht) * inv));
                }
            }
            full[o+3] = (unsigned char)(a2 * 255.0 + 0.5);
        }
        // α=0 の画素の色を、隣の魚の色で埋める。
        // 埋めないまま α をならすと、**背景の白が縁ににじみ出す**
        for (int it = 0; it < 3; ++it) {
            std::vector<unsigned char> cp = full;
            for (int j2 = 1; j2 < SH - 1; ++j2) for (int i2 = 1; i2 < SW - 1; ++i2) {
                size_t o = (size_t)(j2 * SW + i2) * 4;
                if (cp[o+3] > 0) continue;
                int sr = 0, sg = 0, sb = 0, c2 = 0;
                for (int b3 = -1; b3 <= 1; ++b3) for (int a3 = -1; a3 <= 1; ++a3) {
                    size_t o2 = (size_t)((j2 + b3) * SW + (i2 + a3)) * 4;
                    if (cp[o2+3] > 0) { sr += cp[o2+0]; sg += cp[o2+1]; sb += cp[o2+2]; c2++; }
                }
                if (c2) { full[o+0] = (unsigned char)(sr/c2); full[o+1] = (unsigned char)(sg/c2);
                          full[o+2] = (unsigned char)(sb/c2); }
            }
        }
        // 縁を1画素ならす
        std::vector<unsigned char> al2((size_t)SW * SH);
        for (int k = 0; k < SW * SH; ++k) al2[k] = full[(size_t)k * 4 + 3];
        for (int j = 1; j < SH - 1; ++j) for (int i = 1; i < SW - 1; ++i) {
            int sum2 = 0;
            for (int b3 = -1; b3 <= 1; ++b3) for (int a3 = -1; a3 <= 1; ++a3)
                sum2 += al2[(size_t)(j + b3) * SW + (i + a3)];
            full[((size_t)j * SW + i) * 4 + 3] = (unsigned char)(sum2 / 9);
        }
        // 貼るときに要る位置: 鼻先と尾柄がスプライトのどこか、体の中心線の高さ、体長との比
        int xn2 = headRight ? (bestB - 1) : bestA, xp2 = headRight ? bestA : (bestB - 1);
        sprNose = (headRight ? (fx1 - xn2) : (xn2 - fx0)) / (fx1 - fx0);
        sprPed  = (headRight ? (fx1 - xp2) : (xp2 - fx0)) / (fx1 - fx0);
        double ay3 = ((tS[xn2] + bS[xn2]) * 0.5 + (tS[xp2] + bS[xp2]) * 0.5) * 0.5;
        sprAxis = (ay3 - fy0) / (fy1 - fy0);
        sprScale = (fx1 - fx0) / bodyLenGlobal;
        sprAspect = (fy1 - fy0) / (fx1 - fx0);
        fprintf(stderr, "まるごと: 鼻先 %.3f 尾柄 %.3f 中心線 %.3f 幅/体長 %.2f 縦横比 %.2f\n",
                sprNose, sprPed, sprAxis, sprScale, sprAspect);
    }

    if (dbg) stbi_write_png(dbg, SW, SH, 4, full.data(), SW * 4);

    FILE* f = fopen(out, "w");
    if (f) {
        fprintf(f, "// 自動生成 (tools/mktex.cpp)。写真から作った金魚のスプライト。\n");
        fprintf(f, "// 体の輪郭。中心線からの張り出し(体長比)。0=鼻先 1=尾柄\n");
        fprintf(f, "#define %s_HAS_FULL 1\n", sym);
        fprintf(f, "static const int %s_FULL_W = %d, %s_FULL_H = %d;\n", sym, SW, sym, SH);
        fprintf(f, "static const float %s_FULL_MET[5] = {%.4ff,%.4ff,%.4ff,%.4ff,%.4ff};\n",
                sym, sprNose, sprPed, sprAxis, sprScale, sprAspect);
        fprintf(f, "static const unsigned char %s_FULL[] = {\n", sym);
        for (size_t i = 0; i < full.size(); ++i)
            fprintf(f, "%d,%s", full[i], (i % 32 == 31) ? "\n" : "");
        fprintf(f, "};\n");
        fprintf(f, "#define %s_HAS_PROFILE 1\n", sym);
        fprintf(f, "static const int %s_NSEG = %d;\n", sym, NSEG);
        fprintf(f, "static const float %s_MID[] = {", sym);
        for (int i2 = 0; i2 <= NSEG; ++i2) fprintf(f, "%.4ff,", profT[i2]);
        fprintf(f, "};\nstatic const float %s_BOT[] = {", sym);
        for (int i2 = 0; i2 <= NSEG; ++i2) fprintf(f, "%.4ff,", profB[i2]);
        fprintf(f, "};\n");
        fclose(f);
        fprintf(stderr, "wrote %s (%dx%d)\n", out, SW, SH);
    }
    stbi_image_free(img);
    return 0;
}
