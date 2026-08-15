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
static const double INSET = 0.12;   // 背びれ・尻びれの付け根(明るい)を拾わない程度に内側へ

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
        // 尾柄の端は**尾びれの根元(明るい黄色)**が混ざる。そのまま貼ると、体と尾の
        // 継ぎ目に白っぽい塊が出る。最後の1割は 90% の位置の色をそのまま伸ばす
        double uu = std::min(u, 0.90);
        double xf = headRight ? (bestB - 1 - uu * (bestB - 1 - bestA))
                              : (bestA + uu * (bestB - 1 - bestA));
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

    const double bodyLenGlobal = (double)(bestB - bestA);
    // 胴体の平均色。ひれの色みをここに寄せる(ひれと体は同じ色素なので、
    // 白背景ぶんを引いただけの色は黄色に振り切れてしまう)
    double bodyR = 0, bodyG = 0, bodyB = 0;
    {
        long c3 = 0;
        for (int y = 0; y < h; ++y) for (int x = 0; x < w; ++x)
            if (bodym[(size_t)y * w + x]) { bodyR += A(x,y,0); bodyG += A(x,y,1); bodyB += A(x,y,2); c3++; }
        if (c3 < 1) c3 = 1;
        bodyR /= c3; bodyG /= c3; bodyB /= c3;
        fprintf(stderr, "胴体の平均色 (%.0f,%.0f,%.0f)\n", bodyR, bodyG, bodyB);
    }

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

    // ---- ひれも写真から切り出す。**手続きで描いた板では本物に見えない。**
    // 連結成分では分けられない(体の明るいところでひれ同士が繋がる)ので、
    // 胴体の輪郭を基準に「どこにあるか」で分ける:
    //   尾   … 胴体の x 範囲より尾側
    //   背   … 胴体の範囲内で上の縁より外
    //   腹側 … 下の縁より外を、前から 胸/腹/尻 に切る
    // 切り出したら**根元が左、先が右**に揃えて長方形へ詰める(貼る側の向きを1つにするため)。
    // 形は α が持つので、まっすぐな縁の板にはならない。
    struct Fin { const char* name; std::vector<unsigned char> img; int w, h;
                 double lenBL, spanBL, atS; };
    std::vector<Fin> fins;
    const int FW2 = 128, FH2 = 96;
    double bodyLen = bodyLenGlobal;
    auto cut = [&](const char* name, int x0, int x1, bool aboveTop, bool belowBot, int rot) {
        // rot: 0=そのまま(尾) 1=下が根元(背びれ) 2=上が根元(腹側のひれ)
        int ax = w, bx = -1, ay = h, by = -1; long cnt = 0;
        auto isFinRaw = [&](int x, int y) {
            if (x < 0 || y < 0 || x >= w || y >= h) return false;
            if (!fishm[(size_t)y * w + x] || bodym[(size_t)y * w + x]) return false;
            bool inBody = (x >= bestA && x < bestB);
            if (aboveTop && !(inBody && y < tS[x] - 1)) return false;
            if (belowBot && !(inBody && y > bS[x] + 1)) return false;
            return true;
        };
        // **いちばん大きい塊だけ残す。** 写真には髭のような細い線や写り込みが入っていて、
        // そのまま貼ると魚から黒い筋が生えたように見える。
        // ただし細い線はひれと**繋がっている**ことがあるので、塊を数える前に
        // 2画素ぶん痩せさせる(収縮)。細い線はここで消え、ひれは残る
        std::vector<unsigned char> thin((size_t)w * h, 0);
        for (int y = 0; y < h; ++y) for (int x = std::max(0,x0); x <= std::min(w-1,x1); ++x)
            thin[(size_t)y * w + x] = isFinRaw(x, y) ? 1 : 0;
        for (int it2 = 0; it2 < 2; ++it2) {
            std::vector<unsigned char> t2 = thin;
            for (int y = 1; y < h - 1; ++y) for (int x = std::max(1,x0); x <= std::min(w-2,x1); ++x) {
                if (!t2[(size_t)y * w + x]) continue;
                if (!t2[(size_t)(y-1) * w + x] || !t2[(size_t)(y+1) * w + x] ||
                    !t2[(size_t)y * w + x - 1] || !t2[(size_t)y * w + x + 1])
                    thin[(size_t)y * w + x] = 0;
            }
        }
        std::vector<int> lab((size_t)w * h, 0);
        int nlab = 0, bestLab = 0; long bestCnt = 0;
        std::vector<int> st;
        for (int y = 0; y < h; ++y) for (int x = std::max(0,x0); x <= std::min(w-1,x1); ++x) {
            if (!thin[(size_t)y * w + x] || lab[(size_t)y * w + x]) continue;
            ++nlab; long c2 = 0; st.clear(); st.push_back(y * w + x);
            lab[(size_t)y * w + x] = nlab;
            while (!st.empty()) {
                int q2 = st.back(); st.pop_back(); c2++;
                int qx = q2 % w, qy = q2 / w;
                const int dx4[4] = {1,-1,0,0}, dy4[4] = {0,0,1,-1};
                for (int d2 = 0; d2 < 4; ++d2) {
                    int nx2 = qx + dx4[d2], ny2 = qy + dy4[d2];
                    if (nx2 < std::max(0,x0) || nx2 > std::min(w-1,x1) || ny2 < 0 || ny2 >= h) continue;
                    if (lab[(size_t)ny2 * w + nx2] || !thin[(size_t)ny2 * w + nx2]) continue;
                    lab[(size_t)ny2 * w + nx2] = nlab; st.push_back(ny2 * w + nx2);
                }
            }
            // 面積だけで選ぶと、写真に写り込んだ**細い線**が勝ってしまう(長いので画素数が多い)。
            // ひれは面で埋まっているので、外接矩形の詰まり具合を掛けて選ぶ
            long bx0 = w, bx1 = -1, by0 = h, by1 = -1;
            for (int yy = 0; yy < h; ++yy) for (int xx = std::max(0,x0); xx <= std::min(w-1,x1); ++xx)
                if (lab[(size_t)yy * w + xx] == nlab) {
                    bx0 = std::min(bx0, (long)xx); bx1 = std::max(bx1, (long)xx);
                    by0 = std::min(by0, (long)yy); by1 = std::max(by1, (long)yy);
                }
            double fill = (bx1 >= bx0) ? (double)c2 / ((bx1-bx0+1.0) * (by1-by0+1.0)) : 0.0;
            double score = c2 * fill;
            if (score > bestCnt) { bestCnt = (long)score; bestLab = nlab; }
        }
        // 選んだ塊の外接矩形を少し広げ、その中では元のマスクを使う(痩せた分を戻す)
        int kx0 = w, kx1 = -1, ky0 = h, ky1 = -1;
        for (int y = 0; y < h; ++y) for (int x = std::max(0,x0); x <= std::min(w-1,x1); ++x)
            if (lab[(size_t)y * w + x] == bestLab) {
                kx0 = std::min(kx0,x); kx1 = std::max(kx1,x);
                ky0 = std::min(ky0,y); ky1 = std::max(ky1,y);
            }
        const int PAD = 3;
        kx0 -= PAD; kx1 += PAD; ky0 -= PAD; ky1 += PAD;
        auto isFin = [&](int x, int y) {
            if (x < kx0 || x > kx1 || y < ky0 || y > ky1) return false;
            return isFinRaw(x, y);
        };
        for (int y = 0; y < h; ++y) for (int x = std::max(0,x0); x <= std::min(w-1,x1); ++x)
            if (isFin(x, y)) { ax = std::min(ax,x); bx = std::max(bx,x);
                               ay = std::min(ay,y); by = std::max(by,y); cnt++; }
        if (cnt < 60) { fprintf(stderr, "  %-8s 見つからない(%ld画素)\n", name, cnt); return; }
        int cw = bx - ax + 1, ch = by - ay + 1;
        Fin fn; fn.name = name; fn.w = FW2; fn.h = FH2;
        fn.img.assign((size_t)FW2 * FH2 * 4, 0);
        // **根元を揃えて詰める。** 外接矩形をそのまま貼ると、根元の行がほとんど透明な
        // ひれ(基底が斜めのもの)は体から浮いて見える。列ごとに「体に一番近い画素」を探し、
        // そこを u=0 に揃える。長さは列ごとに違ってよい(先の形は α が持つ)
        int maxLen = 1;
        std::vector<int> rootAt(FW2, -1), lenAt(FW2, 0);
        for (int j = 0; j < FH2; ++j) {
            double v = (double)j / (FH2 - 1);
            int x, first = -1, last = -1;
            if (rot == 0) {                       // 尾: 体側は headRight なら大きい x
                int y = (int)(ay + v * (ch - 1) + 0.5);
                y = std::min(h-1, std::max(0, y));
                for (int k = 0; k < cw; ++k) {
                    x = headRight ? (bx - k) : (ax + k);
                    if (isFin(x, y)) { if (first < 0) first = k; last = k; }
                }
            } else {                              // 背/腹: 体側は rot==1 なら大きい y
                x = (int)((headRight ? (bx - v * (cw - 1)) : (ax + v * (cw - 1))) + 0.5);
                x = std::min(w-1, std::max(0, x));
                for (int k = 0; k < ch; ++k) {
                    int y = (rot == 1) ? (by - k) : (ay + k);
                    if (isFin(x, y)) { if (first < 0) first = k; last = k; }
                }
            }
            rootAt[j] = first; lenAt[j] = (first < 0) ? 0 : (last - first + 1);
        }
        if (rot == 0) {
            // **尾は根元が縦一直線**。行ごとに揃えると二叉の切れ込みが潰れてしまう。
            // いちばん体に近い位置を全行の共通の根元にする
            int mn = 1 << 30;
            for (int j = 0; j < FH2; ++j) if (rootAt[j] >= 0) mn = std::min(mn, rootAt[j]);
            if (mn == (1 << 30)) mn = 0;
            for (int j = 0; j < FH2; ++j) {
                if (rootAt[j] < 0) continue;
                lenAt[j] += rootAt[j] - mn;
                rootAt[j] = mn;
            }
        }
        for (int j = 0; j < FH2; ++j) maxLen = std::max(maxLen, lenAt[j]);
        for (int j = 0; j < FH2; ++j) for (int i = 0; i < FW2; ++i) {
            double u = (double)i / (FW2 - 1), v = (double)j / (FH2 - 1);
            size_t o = (size_t)(j * FW2 + i) * 4;
            if (rootAt[j] < 0) { fn.img[o+3] = 0; continue; }
            int step = rootAt[j] + (int)(u * (maxLen - 1) + 0.5);
            int x, y;
            if (rot == 0) {
                x = headRight ? (bx - step) : (ax + step);
                y = (int)(ay + v * (ch - 1) + 0.5);
            } else {
                x = (int)((headRight ? (bx - v * (cw - 1)) : (ax + v * (cw - 1))) + 0.5);
                y = (rot == 1) ? (by - step) : (ay + step);
            }
            x = std::min(w-1, std::max(0, x)); y = std::min(h-1, std::max(0, y));
            if (!isFin(x, y)) { fn.img[o+3] = 0; continue; }
            // **ひれは透ける。** 白背景で撮った写真では、透けたぶんだけ背景の白が
            // 混ざって「白っぽい」色になっている。そのまま貼ると、透けない板になる。
            // 白がどれだけ混ざったかから透過率を出して、ひれ本来の色に戻す:
            //   観測 = a*色 + (1-a)*白   →   a = 1 - min(観測)/255,  色 = (観測 - (1-a)*白)/a
            double r0 = A(x,y,0), g0 = A(x,y,1), b0 = A(x,y,2);
            double mn = std::min(r0, std::min(g0, b0));
            double al = 1.0 - mn / 255.0;
            if (al < 0.06) { fn.img[o+3] = 0; continue; }      // ほぼ透明なところ
            double inv = 1.0 / al, wht = (1.0 - al) * 255.0;
            double cr = (r0 - wht) * inv, cg = (g0 - wht) * inv, cb = (b0 - wht) * inv;
            // **暗い画素は「透けたひれ」ではない。** 体の影が落ちているところや写真の
            // 写り込みで、白背景の前提が崩れている。戻した色が黒くなるので落とす
            // (残すと、ひれの付け根にゴミのような黒い点として出る)
            if (cr + cg + cb < 250.0) { fn.img[o+3] = 0; continue; }
            // **色みは体に合わせる。** 引き算だけだと青が 0 に張り付いて真っ黄色になる。
            // 明暗のむら(鰭条)はひれ自身のものを残し、色みだけ胴体の平均色へ寄せる
            {
                double lum = 0.30 * cr + 0.59 * cg + 0.11 * cb;
                double bl2 = 0.30 * bodyR + 0.59 * bodyG + 0.11 * bodyB;
                double k2 = (bl2 > 1.0) ? lum / bl2 : 1.0;
                const double MIX = 0.62;                  // 体の色みへ寄せる割合
                cr = cr * (1 - MIX) + bodyR * k2 * MIX;
                cg = cg * (1 - MIX) + bodyG * k2 * MIX;
                cb = cb * (1 - MIX) + bodyB * k2 * MIX;
            }
            fn.img[o+0] = (unsigned char)std::min(255.0, std::max(0.0, cr));
            fn.img[o+1] = (unsigned char)std::min(255.0, std::max(0.0, cg));
            fn.img[o+2] = (unsigned char)std::min(255.0, std::max(0.0, cb));
            // 実際のひれはここまで薄くない(白背景ぶんを引きすぎる)ので、少し持ち上げる。
            // 根元は体の下に隠れているので、そこだけ消していく
            double rootFade = std::min(1.0, u / 0.06);
            fn.img[o+3] = (unsigned char)std::min(255.0, (60.0 + al * 210.0) * rootFade);
        }
        // 縁を1画素ならす(切り口のギザギザを消す)
        std::vector<unsigned char> al((size_t)FW2 * FH2);
        for (int k = 0; k < FW2 * FH2; ++k) al[k] = fn.img[(size_t)k * 4 + 3];
        for (int j = 1; j < FH2 - 1; ++j) for (int i = 1; i < FW2 - 1; ++i) {
            int sum = 0;
            for (int b2 = -1; b2 <= 1; ++b2) for (int a2 = -1; a2 <= 1; ++a2)
                sum += al[(size_t)(j + b2) * FW2 + (i + a2)];
            fn.img[((size_t)j * FW2 + i) * 4 + 3] = (unsigned char)(sum / 9);
        }
        fn.lenBL  = maxLen / bodyLen;
        fn.spanBL = (rot == 0 ? ch : cw) / bodyLen;
        double cx2 = (ax + bx) * 0.5;
        fn.atS = headRight ? (bestB - 1 - cx2) / bodyLen : (cx2 - bestA) / bodyLen;
        if (fn.atS < 0) fn.atS = 0; if (fn.atS > 1) fn.atS = 1;
        fins.push_back(fn);
        fprintf(stderr, "  %-8s %3dx%-3d 画素  長さ %.2f体長 幅 %.2f  s=%.2f\n",
                name, cw, ch, fn.lenBL, fn.spanBL, fn.atS);
    };
    fprintf(stderr, "ひれを切り出す:\n");
    cut("CAUDAL",   headRight ? 0 : bestB, headRight ? bestA : w - 1, false, false, 0);
    cut("DORSAL",   bestA, bestB - 1, true, false, 1);
    {
        int L = bestB - bestA;
        int f0 = headRight ? bestB - 1 - (int)(0.42 * L) : bestA + (int)(0.42 * L);
        int f1 = headRight ? bestB - 1 - (int)(0.62 * L) : bestA + (int)(0.62 * L);
        cut("PECTORAL", headRight ? f0 : bestA, headRight ? bestB - 1 : f0, false, true, 2);
        cut("PELVIC",   headRight ? f1 : f0,    headRight ? f0 : f1,        false, true, 2);
        cut("ANAL",     headRight ? bestA : f1, headRight ? f1 : bestB - 1, false, true, 2);
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
        // ひれも写真から。形は α が持つので、まっすぐな縁の板にならない
        for (size_t k = 0; k < fins.size(); ++k) {
            const Fin& fn = fins[k];
            fprintf(f, "static const int %s_%s_W = %d, %s_%s_H = %d;\n",
                    sym, fn.name, fn.w, sym, fn.name, fn.h);
            fprintf(f, "// 長さ/体長, 幅/体長, 付く位置 s(0=鼻先 1=尾柄)\n");
            fprintf(f, "static const float %s_%s_MET[3] = {%.4ff,%.4ff,%.4ff};\n",
                    sym, fn.name, fn.lenBL, fn.spanBL, fn.atS);
            fprintf(f, "static const unsigned char %s_%s[] = {\n", sym, fn.name);
            for (size_t i = 0; i < fn.img.size(); ++i)
                fprintf(f, "%d,%s", fn.img[i], (i % 32 == 31) ? "\n" : "");
            fprintf(f, "};\n");
        }
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
        fprintf(f, "#define %s_NFIN %d\n", sym, (int)fins.size());
        fclose(f);
        fprintf(stderr, "wrote %s (%dx%d, ひれ %d枚)\n", out, TW, TH, (int)fins.size());
    }
    stbi_image_free(img);
    return 0;
}
