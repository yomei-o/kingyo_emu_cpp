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

static const int TW = 256, TH = 80;    // 出力テクスチャ
// 輪郭ぴったりだと背景やひれの根元を拾うので、少し内側に詰めてから展開する
static const double INSET = 0.05;

int main(int argc, char** argv) {
    const char* in   = argc > 1 ? argv[1] : "kingyo.png";
    const char* out  = argc > 2 ? argv[2] : "src/kingyo_tex.h";
    const char* dbg  = argc > 3 ? argv[3] : nullptr;   // 展開結果の確認用
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

    // --- 展開
    std::vector<unsigned char> tex((size_t)TW * TH * 4, 0);
    for (int i = 0; i < TW; ++i) {
        double u = (double)i / (TW - 1);                  // 0=鼻先 1=尾柄
        double xf = headRight ? (bestB - 1 - u * (bestB - 1 - bestA))
                              : (bestA + u * (bestB - 1 - bestA));
        int x = std::min(bestB - 1, std::max(bestA, (int)(xf + 0.5)));
        double t = tS[x], b = bS[x], mid = (t + b) * 0.5, hh = (b - t) * 0.5 * (1.0 - INSET);
        for (int j = 0; j < TH; ++j) {
            double v = (double)j / (TH - 1) * 2.0 - 1.0;   // -1=背 +1=腹
            double yf = mid + v * hh;
            int y = std::min(h - 1, std::max(0, (int)(yf + 0.5)));
            size_t o = (size_t)(j * TW + i) * 4;
            tex[o+0] = (unsigned char)A(x, y, 0);
            tex[o+1] = (unsigned char)A(x, y, 1);
            tex[o+2] = (unsigned char)A(x, y, 2);
            tex[o+3] = 255;
        }
    }
    if (dbg) stbi_write_png(dbg, TW, TH, 4, tex.data(), TW * 4);
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

    FILE* f = fopen(out, "w");
    if (f) {
        fprintf(f, "// 自動生成 (tools/mktex.cpp)。写真から展開した胴体テクスチャ。\n");
        fprintf(f, "static const int %s_TW = %d;\n", sym, TW);
        fprintf(f, "static const int %s_TH = %d;\n", sym, TH);
        // ひれの色も写真から取る(胴体と分けた「明るいほう」の平均)。
        // シミュレータ側はこれを薄めてひれを塗る。手で決めた色だと魚を変えるたびに合わなくなる
        {
            long fr2 = 0, fg2 = 0, fb2 = 0, fc = 0;
            for (int y = 0; y < h; ++y) for (int x = 0; x < w; ++x)
                if (fishm[(size_t)y * w + x] && !bodym[(size_t)y * w + x]) {
                    fr2 += A(x,y,0); fg2 += A(x,y,1); fb2 += A(x,y,2); fc++;
                }
            if (fc < 1) fc = 1;
            fprintf(f, "static const unsigned char %s_FIN[3] = {%d,%d,%d};\n",
                    sym, (int)(fr2/fc), (int)(fg2/fc), (int)(fb2/fc));
            fprintf(stderr, "ひれの平均色 (%d,%d,%d)\n", (int)(fr2/fc), (int)(fg2/fc), (int)(fb2/fc));
        }
        fprintf(f, "static const unsigned char %s_TEX[] = {\n", sym);
        for (size_t i = 0; i < tex.size(); ++i)
            fprintf(f, "%d,%s", tex[i], (i % 32 == 31) ? "\n" : "");
        fprintf(f, "};\n");
        fclose(f);
        fprintf(stderr, "wrote %s (%dx%d)\n", out, TW, TH);
    }
    stbi_image_free(img);
    return 0;
}
