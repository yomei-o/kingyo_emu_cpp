// Kingyo OS — 金魚水槽エミュレータ / 和金   (C++ / WASM)
//
// 魚のアニメーションを「描く」のではなく、泳ぎの物理を積分している。
//
//  1) 遊泳: 魚は体を進行波で送って進む(subcarangiform swimming)。体の中心線の横ずれは
//        h(s,t) = A(s) sin(2π s/λ − 2π f t),   A(s) = A_max s^2.2   (尾に向かって増大)
//     これが「生きている感じ」の正体で、体の形を決めているのはこの1本の式だけ。
//  2) 尾の振り数は速度から決まる。魚類遊泳のストローハル数はどの種でもほぼ一定で
//        St = f·(2A)/U ≈ 0.3   →   f = St·U/(2A)
//     なので速く泳げば勝手に速く振れる。振り数を手で決めてはいない。
//  3) 群れ: 分離・整列・結合(boids)＋壁からの反発＋餌への誘引。
//     向きを変える角速度には上限があり、曲がっている間は体が定常的に湾曲する。
//  4) 体色: ネオンテトラの青緑のラインは色素ではなく構造色(虹色素胞)なので、
//     体側をどの角度から見るかで明るさと色相が変わる。そこも視線との角度から出している。
//  5) 描画: 魚は自前の縦スパン・ラスタライザで塗る。輪郭は被覆率でアンチエイリアスされ、
//     体色は縦方向の連続階調で出る(帯のベタ塗りにしない)。
//     背景は視線と水槽(箱)の内側との交点で面を決める。カメラは固定なので前計算しておき、
//     毎フレームは低解像度で作った caustics と光条を足すだけにしてある。
//
// 計算も描画もすべて C++。olive.c は三角形塗りとテキストだけに使っている。
#define OLIVEC_IMPLEMENTATION
#include "olive.c"
// 胴体は写真から展開したテクスチャを貼る(tools/mktex.cpp が生成)。
// 手元に src/kingyo_tex.h があるときだけ写真を使い、無ければ手続きで作った更紗模様で描く。
#if __has_include("kingyo_tex.h")
  #include "kingyo_tex.h"
  #define HAVE_PHOTO 1
#else
  #define HAVE_PHOTO 0
#endif
#include <vector>
#include <cstdint>
#include <cmath>
#include <cstring>
#include <algorithm>
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define KEEP EMSCRIPTEN_KEEPALIVE
#else
#define KEEP
#endif

// ---------------------------------------------------------------- 画面と水槽
static const int FW = 960, FH = 600;
static const int SUBSTEPS = 2;
static const double DT = 1.0 / (60.0 * SUBSTEPS);

// 水の屈折率。ガラス越しに水中を見ると、ものは実際の 1/n の深さに見える
// (見かけの深さ = 実深さ / n)。水槽が本来より浅く=遠近が弱く見えるのはこれ。
static const double NWATER = 1.333;
// 奥行きは 30cm。魚の前後の重なりと水の減衰はこれで決まる。
static const double TD = 0.30;

// 画角。前面ガラスのところで横 VIEWW[m] が画面いっぱいに写る。
// カメラを引いて焦点距離を伸ばすと、画角は同じまま遠近だけ弱まる(望遠の圧縮)。
// 前後の大きさの比は (TD/n + CAMD) / CAMD なので、0.4m だと 1.56倍、3.0m だと 1.075倍。
static double VIEWW = 0.38;
static const double VIEWW_MAX = 1.50;            // 45インチのテレビで等身大に近い引きまで
static double CAMD = 3.00;                       // 前面ガラスからカメラまで [m]
static const double CAMD_MIN = 0.40;             // スライダの下限(近いほど遠近が強く出る)

// 水槽の幅と高さは決め打ちにせず、**いまの画角から作る**。
// 「床・水面・左右のガラスが画面に入らない」だけが条件で、視線は奥へ行くほど広がるので
// 効くのは一番奥の面での広がりだけ。そこでの倍率は (CAMD + TD/n) / CAMD
// (カメラが近いほど大きい: 0.4m で 1.56倍、3.0m で 1.075倍)。
//
// 条件は「これ以上」なので、そこにさらに WING だけ水を足す。ぴったりに作ると
// **魚が画面の外に出られない**。見えている範囲が水槽そのものだと、フレームの縁で
// 見えない壁に当たって引き返すので、スノードームの中のように見えてしまう。
// 魚(3cm)の 6体長ぶん余らせておけば、画面の外へ抜けて別の個体が入ってくる。
static const double WING = 0.18;                 // 視界の外に残す水 [m]
static double TW = 0.77, TH = 0.62;              // 水槽の幅・高さ [m] (奥行きは TD 固定)
static double CAMY = 0.31;                       // カメラの高さ = 水槽の中央
static double FOC = FW * 3.00 / 0.38;

static void tank_from_view() {
    double mag  = (CAMD + TD / NWATER) / CAMD;
    double halfW = VIEWW * 0.5 * mag;                        // 一番奥での見え幅の半分
    double halfH = VIEWW * ((double)FH / FW) * 0.5 * mag;
    TW   = 2.0 * (halfW + WING);
    TH   = 2.0 * (halfH + WING);
    CAMY = halfH + WING;                                     // 床は視界の下端より WING 下
}

// 画面のこの深さでの見え幅の半分 [m]。餌を「見えているところ」に落とすのに使う
static inline double view_half_w(double z) { return (FW * 0.5) / FOC * (z / NWATER + CAMD); }
static inline double view_half_h(double z) { return (FH * 0.5) / FOC * (z / NWATER + CAMD); }

static const int CW = FW / 4, CH = FH / 4;       // 光の低解像度バッファ
// 水槽の大きさが画角で変わるので、匹数ではなく**密度[匹/L]**を持つ。
// 匹数を固定にすると、寄ると詰まり、引くとがらがらになる。
// 上限は負荷の頭打ち(画角1.5・カメラ0.4m だと 1480L ≒ 1185匹でここに当たる)。
static const size_t MAX_FISH = 1200;

// ---------------------------------------------------------------- 乱数
static uint32_t rng = 22222227u;
static inline double rnd() { rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5; return (rng & 0xFFFFFF) / (double)0x1000000; }
static inline double rnds() { return rnd() * 2.0 - 1.0; }
static inline int clampi(int v, int a, int b) { return v < a ? a : (v > b ? b : v); }
static inline uint32_t rgb(int r, int g, int b) {
    return 0xFF000000u | ((uint32_t)clampi(b,0,255) << 16) | ((uint32_t)clampi(g,0,255) << 8) | (uint32_t)clampi(r,0,255);
}
static inline uint32_t rgba(int r, int g, int b, int a) {
    return ((uint32_t)clampi(a,0,255) << 24) | ((uint32_t)clampi(b,0,255) << 16)
         | ((uint32_t)clampi(g,0,255) << 8) | (uint32_t)clampi(r,0,255);
}
static inline double smoothstep(double a, double b, double x) {
    if (b <= a) return x < a ? 0.0 : 1.0;
    double t = (x - a) / (b - a);
    if (t < 0) t = 0; if (t > 1) t = 1;
    return t * t * (3.0 - 2.0 * t);
}

// ---------------------------------------------------------------- 投影
static inline void project(double x, double y, double z, double& sx, double& sy, double& sc) {
    double zz = z / NWATER + CAMD;          // 屈折で奥行きが 1/n に縮んで見える
    if (zz < 0.05) zz = 0.05;
    sc = FOC / zz;
    sx = FW * 0.5 + x * sc;
    sy = FH * 0.5 - (y - CAMY) * sc;
}

// ---------------------------------------------------------------- 魚
struct Fish {
    double x, y, z;        // 位置(体の中心) [m]
    double yaw, pitch;     // 向き [rad]
    double speed;          // 遊泳速度 [m/s]
    double L;              // 全長 [m]
    double phase;          // 尾を振る位相 [rad]
    double freq;           // 振り数 [Hz] — 速度から決まる
    double bend;           // 旋回による体の定常的な湾曲
    double dartT;          // 突進(スパート)の残り時間 [s]
    // burst-and-coast: 数回振って加速 → 体をまっすぐにして滑空、の繰り返し
    int    gait;           // 0=バースト(振っている) 1=滑空(止めている)
    double beats;          // このバーストで残り何回振るか
    double coastT;         // 滑空の残り時間 [s]
    double ampS;           // 尾の振幅スケール 0..1 (滑空では 0 に落ちる)
    double cruiseK;        // 個体ごとの巡航速度のくせ
    double girth;          // 体高の倍率(痩せている個体、太っている個体)
    float  tint;           // 個体差(明るさ)
    float  pat;            // 更紗(赤白斑)の位置。個体ごとにずらす
    float  sara;           // 0=素赤(全身朱) 〜 1=白勝ち更紗
    double pecPh;          // 胸びれを漕ぐ位相
};
static std::vector<Fish> fish;

struct Food { double x, y, z, t; };
static std::vector<Food> foods;
static double startle = 0.0;

static std::vector<uint32_t> px;       // 出力 RGBA
static std::vector<uint32_t> bgStatic; // 動かない背景(カメラ固定なので前計算できる)
static std::vector<float>    caU, caV, caW;   // 低解像度セルの世界座標
static std::vector<uint8_t>  caFace;
static std::vector<float>    lightBuf;        // 低解像度の追加光 (RGB)
static bool bgAllBack = false;                // 全セルが奥のガラス = 光を列と行に分けられる

// ---------------------------------------------------------------- パラメータ
static double p_dens   = 0.05;       // 匹/L (金魚は 1匹あたり 10〜20L。20Lに1匹で 0.05)
static double p_school = 1.0;
static double p_speed  = 1.0;
static double p_murk   = 1.0;
static double p_light  = 1.0;
static double p_hud    = 1.0;
// モデル確認モード: 1匹だけを画面中央に固定して泳がせる(体は波打ち、ひれもなびく)。
// 形やひれを直すときに、魚が画面に入るのを待たなくて済む
static bool   p_model  = false;

static double simTime = 0.0;
static long   frameNo = 0;

// ---------------------------------------------------------------- 体の形
// 全長に対する体高の比を表で持つ。**和金は体高/体長 ≒ 0.36** でテトラ(0.27)よりずっと深い。
// 頭から肩にかけて一気に立ち上がり、背びれの前(3〜4割)で最大、尾柄は太いまま終わる
// (金魚は尾柄が細くくびれない。ここを細くするとテトラの体型に戻ってしまう)。
static const double HPROF[] = { 0.10, 0.52, 0.78, 0.93, 1.00, 1.00, 0.96, 0.89, 0.78, 0.65, 0.50, 0.36, 0.24 };
static const int    NHP = (int)(sizeof(HPROF) / sizeof(HPROF[0]));
static inline double body_height(double s) {          // s: 0(鼻先)〜1(尾柄)
    double u = s * (NHP - 1);
    int i = (int)u; if (i < 0) i = 0; if (i > NHP - 2) i = NHP - 2;
    double f = u - i;
    f = f * f * (3.0 - 2.0 * f);                      // なめらかに繋ぐ
    return HPROF[i] * (1 - f) + HPROF[i + 1] * f;
}
// 体高比と体の厚み(横から見たときの「幅」)。金魚は寸胴で厚みがある
static const double BODY_H = 0.36;    // 体高 / 体長
static const double BODY_W = 0.62;    // 体の厚み / 体高

// 体の中心線を伝わる進行波。**s は 1 を超えて尾びれまで続く。**
// ひれは体の延長で、波はそのまま尾びれへ伝わっていく。s が大きいほど位相が遅れるので、
// 尾びれは体に**遅れてなびく**。金魚らしさはほとんどこれで決まる
// (剛体の扇を尾柄に貼ると、板を振っているようにしか見えない)。
// ひれの中を伝わる波。**ひれは体より柔らかい**ので、同じ長さでも位相が余計に遅れ、
// 先へ行くほど大きく振れる。金魚の尾が「波打つように」曲がるのはこれ。
// q: 0=付け根 1=先端 / lagK: 遅れの強さ / ampK: 先端での振幅の伸び
static inline double fin_wave(double q, double phase, double amp, double lagK, double ampK) {
    double a = amp * (1.0 + ampK * q * q);
    return a * sin(2.0 * M_PI * (1.0 + q * lagK) / 1.15 - phase);
}
static inline double body_wave(double s, double phase, double amp, double bend) {
    double a = amp * pow(s, 2.0);
    return a * sin(2.0 * M_PI * s / 1.15 - phase) + bend * s * s;
}

// ---------------------------------------------------------------- caustics
static inline double caustic(double u, double v, double t) {
    double c = sin(u * 21.0 + t * 1.7)
             + sin(v * 17.0 - t * 1.3 + u * 8.0)
             + sin((u + v) * 13.0 + t * 2.1)
             + sin((u - v) * 27.0 - t * 0.9);
    double b = c * 0.125 + 0.5;
    double b2 = b * b, b4 = b2 * b2;
    return b4 * b;                                     // b^5 (pow より速い)
}

// 視線と箱の内側の交点。face: 0=奥 1=左 2=右 3=床 4=水面
// 水槽は「画角の上限でも床・水面・左右が入らない」大きさにしてあるので、
// ふつうは 0(奥のガラス)しか返らない。残りは、上限を超える画角を直に渡したときの保険。
static inline void ray_box(double dx, double dy, double& wx, double& wy, double& wz, int& face) {
    // 見かけの奥行き(= 実奥行き/n)の空間で交点を出し、最後に実奥行きへ戻す
    double tBest = TD / NWATER + CAMD; face = 0;
    if (dx < -1e-9) { double t = (-TW * 0.5) / dx; if (t > 0 && t < tBest) { tBest = t; face = 1; } }
    if (dx >  1e-9) { double t = ( TW * 0.5) / dx; if (t > 0 && t < tBest) { tBest = t; face = 2; } }
    if (dy < -1e-9) { double t = (0.0 - CAMY) / dy; if (t > 0 && t < tBest) { tBest = t; face = 3; } }
    if (dy >  1e-9) { double t = (TH  - CAMY) / dy; if (t > 0 && t < tBest) { tBest = t; face = 4; } }
    wx = tBest * dx; wy = CAMY + tBest * dy; wz = (tBest - CAMD) * NWATER;
}

// カメラが動かないので、面の判定・砂利・奥行き減衰は一度きりで済む
static void build_static_bg() {
    bgStatic.assign((size_t)FW * FH, 0);
    caU.assign(CW * CH, 0); caV.assign(CW * CH, 0); caW.assign(CW * CH, 0);
    caFace.assign(CW * CH, 0);
    lightBuf.assign((size_t)CW * CH * 3, 0.0f);

    for (int j = 0; j < FH; ++j) {
        double dy = -((double)j - FH * 0.5) / FOC;
        for (int i = 0; i < FW; ++i) {
            double dx = ((double)i - FW * 0.5) / FOC;
            double wx, wy, wz; int face;
            ray_box(dx, dy, wx, wy, wz, face);
            double dist = (face == 3 || face == 4) ? (wz / TD) : ((wz + 0.02) / TD);
            double dep = 1.0 - 0.60 * std::min(1.0, std::max(0.0, dist)) * p_murk;
            int r, g, b;
            if (face == 4) {                                  // 水面(下から)
                r = (int)(22 * dep); g = (int)(52 * dep); b = (int)(56 * dep);
            } else {                                          // 奥・左右・床
                double vg = 1.0 - 0.5 * (1.0 - wy / TH);
                double sidedim = (face == 0) ? 1.0 : 0.66;
                r = (int)((12 + 18 * vg) * dep * sidedim);
                g = (int)((32 + 44 * vg) * dep * sidedim);
                b = (int)((36 + 48 * vg) * dep * sidedim);
            }
            bgStatic[(size_t)j * FW + i] = rgb(r, g, b);
        }
    }
    // 低解像度セルの世界座標も一度きり
    for (int j = 0; j < CH; ++j) {
        double dy = -(((j + 0.5) * 4.0) - FH * 0.5) / FOC;
        for (int i = 0; i < CW; ++i) {
            double dx = (((i + 0.5) * 4.0) - FW * 0.5) / FOC;
            double wx, wy, wz; int face;
            ray_box(dx, dy, wx, wy, wz, face);
            int k = j * CW + i;
            caU[k] = (float)wx; caV[k] = (float)wy; caW[k] = (float)wz; caFace[k] = (uint8_t)face;
        }
    }
    // 全セルが奥のガラスなら、光は列と行に分けて作れる(draw_water の速い道)
    bgAllBack = true;
    for (int k = 0; k < CW * CH; ++k) if (caFace[k] != 0) { bgAllBack = false; break; }
}

// 光を列と行に分けて作る速い道。
// 水槽を画角から作っているので映るのは奥のガラスだけ。すると **wz は一定で、
// wx は列だけ・wy は行だけ**で決まる。caustics の各項は
//   sin(u·a + v·b + c) = sin(u·a)·cos(v·b+c) + cos(u·a)·sin(v·b+c)
// と分けられるので、列ごと・行ごとに sin/cos を作っておけばセルでは掛け算だけになる。
// 光条は wx だけの関数なので、なおさら列ごとに1回で済む。
// sin の回数: 36000セル×7 → (240列 + 150行)×数回。出る絵は同じ。
static void light_cells_back(double t, double la) {
    const double t8 = t * 0.8;
    const double wz0 = caW[0];                       // 奥のガラス(全セル共通)
    static double c1[CW], cPs[CW], cPc[CW], cRs[CW], cRc[CW], cTs[CW], cTc[CW], cSh[CW];
    for (int i = 0; i < CW; ++i) {
        double u = caU[i] + wz0 * 0.6;
        c1[i] = sin(u * 21.0 + t8 * 1.7);
        double P = u * 8.0, R = u * 13.0, T = u * 27.0;
        cPs[i] = sin(P); cPc[i] = cos(P);
        cRs[i] = sin(R); cRc[i] = cos(R);
        cTs[i] = sin(T); cTc[i] = cos(T);
        double lx = caU[i] * 22.0 + wz0 * 2.0;
        double sh = sin(lx * 1.0 + t * 0.30) + sin(lx * 2.37 - t * 0.21) + sin(lx * 0.61 + t * 0.13);
        sh = std::max(0.0, sh * 0.166 + 0.5);
        double s2 = sh * sh, s4 = s2 * s2;
        cSh[i] = s4 * s4 * 26.0 * p_light;
    }
    static double rQs[CH], rQc[CH], rSs[CH], rSc[CH], rWs[CH], rWc[CH], rSh[CH];
    for (int j = 0; j < CH; ++j) {
        double wy = caV[(size_t)j * CW], v = wy + 0.7;
        double Q = v * 17.0 - t8 * 1.3, S = v * 13.0 + t8 * 2.1, W = -v * 27.0 - t8 * 0.9;
        rQs[j] = sin(Q); rQc[j] = cos(Q);
        rSs[j] = sin(S); rSc[j] = cos(S);
        rWs[j] = sin(W); rWc[j] = cos(W);
        rSh[j] = 0.15 + 0.85 * (wy / TH);
    }
    const double ca = la * 0.45;
    for (int j = 0; j < CH; ++j) {
        double Qs = rQs[j], Qc = rQc[j], Ss = rSs[j], Sc = rSc[j], Ws = rWs[j], Wc = rWc[j];
        double shr = rSh[j];
        float* out = &lightBuf[(size_t)j * CW * 3];
        for (int i = 0; i < CW; ++i) {
            double c = c1[i] + (cPs[i] * Qc + cPc[i] * Qs)
                             + (cRs[i] * Sc + cRc[i] * Ss)
                             + (cTs[i] * Wc + cTc[i] * Ws);
            double b = c * 0.125 + 0.5;
            double b2 = b * b, b4 = b2 * b2;
            double cc = b4 * b * ca;                 // caustic() と同じ b^5
            double sh = cSh[i] * shr;
            out[i * 3    ] = (float)(105 * cc + sh * 0.62);
            out[i * 3 + 1] = (float)(145 * cc + sh);
            out[i * 3 + 2] = (float)(136 * cc + sh * 0.9);
        }
    }
}

// 床や水面も映る場合の一般の道(画角の上限を超える値を直に渡したときの保険)
static void light_cells_general(double t, double la) {
    for (int k = 0; k < CW * CH; ++k) {
        double wx = caU[k], wy = caV[k], wz = caW[k];
        int face = caFace[k];
        double ar = 0, ag = 0, ab = 0;
        if (face == 4) {
            double rp = sin(wx * 34.0 + t * 1.9) + sin(wz * 27.0 - t * 1.4) + sin((wx + wz) * 44.0 + t * 2.6);
            double rip = std::max(0.0, rp * 0.166 + 0.5);
            double r2 = rip * rip;
            double c = (0.18 + 0.82 * r2 * r2) * la;
            ar = 52 * c; ag = 96 * c; ab = 100 * c;
        } else {
            double c = caustic(wx + wz * 0.6, wy + 0.7, t * 0.8) * la * 0.45;
            ar = 105 * c; ag = 145 * c; ab = 136 * c;
        }
        // 光条(水面から差し込む筋)
        double lx = wx * 22.0 + wz * 2.0;
        double sh = sin(lx * 1.0 + t * 0.30) + sin(lx * 2.37 - t * 0.21) + sin(lx * 0.61 + t * 0.13);
        sh = std::max(0.0, sh * 0.166 + 0.5);
        double s2 = sh * sh, s4 = s2 * s2;
        sh = s4 * s4 * (0.15 + 0.85 * (wy / TH)) * 26.0 * p_light;
        ar += sh * 0.62; ag += sh; ab += sh * 0.9;
        lightBuf[k * 3] = (float)ar; lightBuf[k * 3 + 1] = (float)ag; lightBuf[k * 3 + 2] = (float)ab;
    }
}

// 低解像度(1/4)の光を画面へ拡大するときの、行と列の重み。画面座標だけで決まる
static int   upI0[FW], upI1[FW], upJ0[FH], upJ1[FH];
static float upTx[FW], upTy[FH];
static void build_upsample_tables() {
    for (int i = 0; i < FW; ++i) {
        double fx = (double)i / 4.0 - 0.5;
        int i0 = clampi((int)floor(fx), 0, CW - 1);
        double tx = fx - i0; if (tx < 0) tx = 0; if (tx > 1) tx = 1;
        upI0[i] = i0; upI1[i] = clampi(i0 + 1, 0, CW - 1); upTx[i] = (float)tx;
    }
    for (int j = 0; j < FH; ++j) {
        double fy = (double)j / 4.0 - 0.5;
        int j0 = clampi((int)floor(fy), 0, CH - 1);
        double ty = fy - j0; if (ty < 0) ty = 0; if (ty > 1) ty = 1;
        upJ0[j] = j0; upJ1[j] = clampi(j0 + 1, 0, CH - 1); upTy[j] = (float)ty;
    }
}

// 毎フレーム: 低解像度で光を作って、双一次補間で背景に足す
static void draw_water() {
    double t = simTime;
    double la = 0.9 * p_light;
    if (bgAllBack) light_cells_back(t, la); else light_cells_general(t, la);
    // 拡大の重みは画面座標だけで決まる = 毎フレーム作り直す必要がない(upsample_tables)。
    // 縦の補間は行ごとに1本の中間バッファへまとめてしまえば、画素ごとに残るのは
    // 横の補間だけ(4本の掛け算 → 2本)。ここが一番回るところなので効く。
    static float rowL[CW * 3];
    for (int j = 0; j < FH; ++j) {
        const float* A = &lightBuf[(size_t)upJ0[j] * CW * 3];
        const float* B = &lightBuf[(size_t)upJ1[j] * CW * 3];
        const float ty = upTy[j];
        for (int k = 0; k < CW * 3; ++k) rowL[k] = A[k] + (B[k] - A[k]) * ty;
        const uint32_t* src = &bgStatic[(size_t)j * FW];
        uint32_t* dst = &px[(size_t)j * FW];
        for (int i = 0; i < FW; ++i) {
            const float* a = &rowL[upI0[i] * 3];
            const float* b = &rowL[upI1[i] * 3];
            const float tx = upTx[i];
            uint32_t s = src[i];
            int r = (int)(s & 255)        + (int)(a[0] + (b[0] - a[0]) * tx);
            int g = (int)((s >> 8) & 255) + (int)(a[1] + (b[1] - a[1]) * tx);
            int bb= (int)((s >> 16) & 255)+ (int)(a[2] + (b[2] - a[2]) * tx);
            dst[i] = rgb(r, g, bb);
        }
    }
}

// ---------------------------------------------------------------- 描画の下請け
static Olivec_Canvas OC;
static inline void blend_px(int x, int y, int r, int g, int b, double a) {
    if (a <= 0.004 || x < 0 || y < 0 || x >= FW || y >= FH) return;
    if (a > 1) a = 1;
    uint32_t& d = px[(size_t)y * FW + x];
    int dr = d & 255, dg = (d >> 8) & 255, db = (d >> 16) & 255;
    d = rgb((int)(dr + (r - dr) * a), (int)(dg + (g - dg) * a), (int)(db + (b - db) * a));
}
static inline void tri(double x1, double y1, double x2, double y2, double x3, double y3, uint32_t c) {
    olivec_triangle(OC, (int)llround(x1), (int)llround(y1), (int)llround(x2), (int)llround(y2),
                    (int)llround(x3), (int)llround(y3), c);
}

// ---------------------------------------------------------------- ひれを貼る
// 写真から切り出したひれ(α 付き)を、曲がった四角形に貼る。
// 形は α が持つので、まっすぐな縁の板にならない。u=根元0→先1、v=帯の片側0→反対1。
static void tex_tri(double x0,double y0,double u0,double v0,
                    double x1,double y1,double u1,double v1,
                    double x2,double y2,double u2,double v2,
                    const unsigned char* T, int tw, int th,
                    double aK, double shade, double fog) {
    double minx = std::min(x0, std::min(x1, x2)), maxx = std::max(x0, std::max(x1, x2));
    double miny = std::min(y0, std::min(y1, y2)), maxy = std::max(y0, std::max(y1, y2));
    int xa = (int)floor(minx), xb = (int)ceil(maxx), ya = (int)floor(miny), yb = (int)ceil(maxy);
    if (xb < 0 || yb < 0 || xa >= FW || ya >= FH) return;
    xa = std::max(0, xa); ya = std::max(0, ya); xb = std::min(FW - 1, xb); yb = std::min(FH - 1, yb);
    double d = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2);
    if (fabs(d) < 1e-9) return;
    for (int y = ya; y <= yb; ++y) for (int x = xa; x <= xb; ++x) {
        double px2 = x + 0.5, py2 = y + 0.5;
        double l0 = ((y1 - y2) * (px2 - x2) + (x2 - x1) * (py2 - y2)) / d;
        double l1 = ((y2 - y0) * (px2 - x2) + (x0 - x2) * (py2 - y2)) / d;
        double l2 = 1.0 - l0 - l1;
        if (l0 < -0.002 || l1 < -0.002 || l2 < -0.002) continue;
        double u = l0 * u0 + l1 * u1 + l2 * u2, v = l0 * v0 + l1 * v1 + l2 * v2;
        double tu = u * (tw - 1), tv = v * (th - 1);
        int i0 = clampi((int)tu, 0, tw - 2), j0 = clampi((int)tv, 0, th - 2);
        double ft = tu - i0, gt = tv - j0;
        const unsigned char* a = &T[((size_t)j0 * tw + i0) * 4];
        const unsigned char* b = a + 4;
        const unsigned char* c = &T[((size_t)(j0 + 1) * tw + i0) * 4];
        const unsigned char* e = c + 4;
        double w00 = (1-ft)*(1-gt), w10 = ft*(1-gt), w01 = (1-ft)*gt, w11 = ft*gt;
        double al = (a[3]*w00 + b[3]*w10 + c[3]*w01 + e[3]*w11) / 255.0 * aK;
        if (al < 0.01) continue;
        double r = (a[0]*w00 + b[0]*w10 + c[0]*w01 + e[0]*w11) * shade;
        double g = (a[1]*w00 + b[1]*w10 + c[1]*w01 + e[1]*w11) * shade;
        double bl= (a[2]*w00 + b[2]*w10 + c[2]*w01 + e[2]*w11) * shade;
        blend_px(x, y, (int)(r * (1 - fog) + 14 * fog),
                       (int)(g * (1 - fog) + 44 * fog),
                       (int)(bl * (1 - fog) + 48 * fog), al);
    }
}
static void tex_quad(double ax,double ay,double bx,double by,double cx,double cy,double dx2,double dy2,
                     double u0,double u1,double v0,double v1,
                     const unsigned char* T,int tw,int th,double aK,double shade,double fog) {
    tex_tri(ax,ay,u0,v0, bx,by,u1,v0, cx,cy,u1,v1, T,tw,th,aK,shade,fog);
    tex_tri(ax,ay,u0,v0, cx,cy,u1,v1, dx2,dy2,u0,v1, T,tw,th,aK,shade,fog);
}

// ---------------------------------------------------------------- 1匹描く
// 体は「縦スパン + 被覆率」で自前に塗る。上下の輪郭がアンチエイリアスされ、
// 体色も v(背=-1 〜 腹=+1) の連続関数として出せるので、帯のベタ塗りにならない。

#if HAVE_PHOTO
// ---------------------------------------------------------------- 写真1枚を曲げて描く
// 体を短冊に展開し、ひれを別々に貼る作りは、写真の「魚らしさ」を分解してしまう。
// 頭の丸み・ひれの付き方・目の位置がそれぞれ別に近似され、どこかが必ず嘘になる。
// ここでは**写真1枚をそのまま格子に貼り、格子を波で曲げる**。
// 出てくる絵は写真そのもので、曲がり方だけがこちらの計算。
//
// 格子の各点は「スプライト上の位置(u,v)」から:
//   u → s(鼻先0 尾柄1、尾びれは1超) → 体軸に沿った位置と、波による横のずれ
//   v → 体の中心線からの上下のずれ
// さらに**写真から測った輪郭より外にある点はひれ**とみなし、余分に揺らす。
// こうすると胸びれも腹びれも尻びれも、切り出さずに動く。
// 格子の三角形を1枚塗る。**明るさと透け具合も頂点で持って補間する。**
// 面ごとに1つの値にすると、格子の目がそのまま縞になって出る
// (胴の丸みも、ひれの透け方も、格子より細かく変わるので面単位では足りない)
static void warp_tri(const double* X, const double* Y, const double* U, const double* V,
                     const double* S, const double* A, double fog) {
    double minx = std::min(X[0], std::min(X[1], X[2])), maxx = std::max(X[0], std::max(X[1], X[2]));
    double miny = std::min(Y[0], std::min(Y[1], Y[2])), maxy = std::max(Y[0], std::max(Y[1], Y[2]));
    int xa = (int)floor(minx), xb = (int)ceil(maxx), ya = (int)floor(miny), yb = (int)ceil(maxy);
    if (xb < 0 || yb < 0 || xa >= FW || ya >= FH) return;
    xa = std::max(0, xa); ya = std::max(0, ya); xb = std::min(FW - 1, xb); yb = std::min(FH - 1, yb);
    double d = (Y[1] - Y[2]) * (X[0] - X[2]) + (X[2] - X[1]) * (Y[0] - Y[2]);
    if (fabs(d) < 1e-9) return;
    const int tw = KINGYO_FULL_W, th = KINGYO_FULL_H;
    for (int y = ya; y <= yb; ++y) for (int x = xa; x <= xb; ++x) {
        double px2 = x + 0.5, py2 = y + 0.5;
        double l0 = ((Y[1] - Y[2]) * (px2 - X[2]) + (X[2] - X[1]) * (py2 - Y[2])) / d;
        double l1 = ((Y[2] - Y[0]) * (px2 - X[2]) + (X[0] - X[2]) * (py2 - Y[2])) / d;
        double l2 = 1.0 - l0 - l1;
        if (l0 < 0.0 || l1 < 0.0 || l2 < 0.0) continue;
        double u = l0 * U[0] + l1 * U[1] + l2 * U[2], v = l0 * V[0] + l1 * V[1] + l2 * V[2];
        double sh = l0 * S[0] + l1 * S[1] + l2 * S[2], ak = l0 * A[0] + l1 * A[1] + l2 * A[2];
        double tu = u * (tw - 1), tv = v * (th - 1);
        int i0 = clampi((int)tu, 0, tw - 2), j0 = clampi((int)tv, 0, th - 2);
        double ft = tu - i0, gt = tv - j0;
        const unsigned char* a = &KINGYO_FULL[((size_t)j0 * tw + i0) * 4];
        const unsigned char* b = a + 4;
        const unsigned char* c = &KINGYO_FULL[((size_t)(j0 + 1) * tw + i0) * 4];
        const unsigned char* e = c + 4;
        double w00 = (1-ft)*(1-gt), w10 = ft*(1-gt), w01 = (1-ft)*gt, w11 = ft*gt;
        double al = (a[3]*w00 + b[3]*w10 + c[3]*w01 + e[3]*w11) / 255.0 * ak;
        if (al < 0.008) continue;
        double r = (a[0]*w00 + b[0]*w10 + c[0]*w01 + e[0]*w11) * sh;
        double g = (a[1]*w00 + b[1]*w10 + c[1]*w01 + e[1]*w11) * sh;
        double bl= (a[2]*w00 + b[2]*w10 + c[2]*w01 + e[2]*w11) * sh;
        blend_px(x, y, (int)(r * (1 - fog) + 14 * fog),
                       (int)(g * (1 - fog) + 44 * fog),
                       (int)(bl * (1 - fog) + 48 * fog), al);
    }
}

static void draw_fish_full(const Fish& f) {
    const float* M = KINGYO_FULL_MET;      // 鼻先u, 尾柄u, 中心線v, 幅/体長, 縦横比
    const double U0 = M[0], U1 = M[1], VA = M[2];
    const double HB = M[3] * M[4];         // スプライトの高さ / 体長
    const double SEND = 1.0 + (1.0 - U1) / (U1 - U0);   // 画像の右端に当たる s

    double cy = cos(f.yaw), sy = sin(f.yaw);
    double fx = cy * cos(f.pitch), fy = sin(f.pitch), fz = sy * cos(f.pitch);
    double ax = -sy, az = cy;

    double csx, csy, sc0;
    project(f.x, f.y, f.z, csx, csy, sc0);
    double lenpx = f.L * sc0;
    const int NXG = clampi((int)(lenpx / 8.0), 8, 40);   // 画面で小さい魚は粗くてよい
    const int NYG = clampi((int)(lenpx / 15.0), 4, 20);
    static const int MXG = 41, MYG = 21;
    double PX2[MXG][MYG], PY2[MXG][MYG];
    double RND[MXG][MYG], ALP[MXG][MYG];   // 断面の丸みによる陰影 / ひれの透け具合
    double WXc[MXG], WZc[MXG], CSH[MXG];   // 列ごとの世界座標と、体側の向きによる明るさ

    double amp  = 0.070 * f.L * f.ampS;
    double bodyW = 0.9 * f.L;              // 鼻先〜尾柄の長さ [m]
    double dep = std::max(0.0, std::min(1.0, f.z / TD));
    double fog = 0.62 * dep * p_murk;

    for (int i = 0; i <= NXG; ++i) {
        double u = (double)i / NXG;
        double sq = (u - U0) / (U1 - U0);
        double along0 = f.L * (0.46 - sq * 0.9);
        // 体の波 → 尾びれのしなりへ。ひれは体より柔らかいので余計に遅れる
        double q = std::max(0.0, (sq - 1.0) / std::max(1e-6, SEND - 1.0));
        double lat0 = body_wave(std::min(1.0, std::max(0.0, sq)) * 0.9, f.phase, amp, f.bend) * (1.0 - q)
                    + fin_wave(q, f.phase, amp, 1.6, 1.8) * q;
        // 写真から測った輪郭。ここより外はひれ
        double u2 = std::min(1.0, std::max(0.0, sq)) * KINGYO_NSEG;
        int i2 = clampi((int)u2, 0, KINGYO_NSEG - 1);
        double fr3 = u2 - i2;
        double midv  = (KINGYO_MID[i2]  + (KINGYO_MID[i2+1]  - KINGYO_MID[i2])  * fr3) / HB;
        double halfv = (KINGYO_BOT[i2] + (KINGYO_BOT[i2+1] - KINGYO_BOT[i2]) * fr3) / HB;
        if (sq > 1.0) halfv *= std::max(0.15, 1.0 - (sq - 1.0) * 2.0);   // 尾柄から先は体でない
        double swing = sin(f.pecPh + sq * 2.6);        // ひれごとに位相をずらす
        // 体は板ではなく**断面が楕円の筒**。見えているのは手前半分なので、
        // 各点をカメラ側へ膨らませる。こうすると斜めを向いたときに
        // 胴が短くなるだけでなく**幅を持ったまま回る**(平面のままだと紙が回るだけ)
        double wx0 = f.x + fx * along0 + ax * lat0;
        double wz0 = f.z + fz * along0 + az * lat0;
        double toCam = ax * (0.0 - wx0) + az * (-CAMD * NWATER - wz0);
        double bulgeSgn = (toCam >= 0.0) ? 1.0 : -1.0;
        double halfThick = 0.62 * halfv * HB * bodyW * f.girth;   // 体の厚みの半分 [m]
        for (int j = 0; j <= NYG; ++j) {
            double v = (double)j / NYG;
            double dv = v - (VA + midv);
            double out = fabs(dv) - halfv;             // >0 ならひれ
            double along = along0, lat = lat0, dvv = dv;
            if (out > 0.0) {
                double ow = out * HB * bodyW;          // 世界座標でのひれの張り出し [m]
                double k = std::min(1.0, out / 0.22);
                dvv = (dv < 0 ? -1.0 : 1.0) * (halfv + out * (1.0 + 0.13 * swing));  // 扇が開閉する
                along -= 0.28 * ow * k * swing;        // 後ろへ流れる
                lat   += 0.85 * ow * k * swing;        // 横に煽ぐ
            }
            double yoff = -dvv * HB * bodyW * f.girth;
            double t2 = dv / std::max(1e-6, halfv);
            double nz0 = (out > 0.0) ? 0.0 : sqrt(std::max(0.0, 1.0 - t2 * t2));
            double blg = bulgeSgn * halfThick * nz0;   // カメラ側への膨らみ
            double sxp, syp, scp;
            project(f.x + fx * along + ax * (lat + blg), f.y + fy * along + yoff,
                    f.z + fz * along + az * (lat + blg), sxp, syp, scp);
            PX2[i][j] = sxp; PY2[i][j] = syp;
            // 断面は楕円。**写真をそのまま貼ると板に見える**ので、
            // 中心線から縁へ行くほど法線が寝る = 暗くなる ぶんを掛ける。
            // 少し背中寄りに濡れた体の照り返しを足す
            if (out > 0.0) {
                RND[i][j] = 0.80;                      // ひれは平らな膜
                // ひれは薄い膜。**黄色いひれは彩度が高いので、写真の明るさから
                // 透け具合を測ると「不透明」と出てしまう**。張り出しの量から決める
                ALP[i][j] = std::max(0.30, 0.90 - 0.55 * std::min(1.0, out / 0.25));
            } else {
                // 法線は「横向き成分 nz0 + 上下成分 -t2」。上からの光で背が明るく、
                // 縁に行くほど法線が寝て暗くなる。濡れた体の照りは背中寄りに出る
                RND[i][j] = 0.62 + 0.30 * nz0 + 0.12 * std::max(0.0, -t2)
                          + 0.16 * exp(-((t2 + 0.40) * (t2 + 0.40)) / 0.05) * nz0;
                ALP[i][j] = 1.0;
            }
        }
        WXc[i] = f.x + fx * along0 + ax * lat0;
        WZc[i] = f.z + fz * along0 + az * lat0;
    }
    // 列ごとの体側の向き。**曲がった胴は、こちらを向いた側が明るくなる**。
    // 1枚の絵をただ曲げるだけでは出ない「厚みのある体が泳いでいる」感じは、ほぼこれで出る
    for (int i = 0; i <= NXG; ++i) {
        int a2 = std::max(0, i - 1), b2 = std::min(NXG, i + 1);
        double tx = WXc[b2] - WXc[a2], tz = WZc[b2] - WZc[a2];
        double tl = sqrt(tx * tx + tz * tz); if (tl < 1e-9) { tx = fx; tz = fz; tl = 1; }
        double nx = -tz / tl, nz2 = tx / tl;             // 体側の法線(水平面内)
        double ex = 0.0 - WXc[i], ez = -CAMD * NWATER - WZc[i];
        double el = sqrt(ex * ex + ez * ez); if (el < 1e-9) el = 1;
        double sd = fabs((nx * ex + nz2 * ez) / el);     // 1=真横を向いている
        CSH[i] = (0.62 + 0.44 * pow(sd, 0.7)) * f.tint;
    }
    // 正面を向いた魚が線に潰れないようにする(魚は板ではなく厚みがある)
    {
        double dmin = PX2[0][0], dmax = PX2[0][0], cx = 0.0; int cnt = 0;
        for (int i = 0; i <= NXG; ++i) for (int j = 0; j <= NYG; ++j) {
            dmin = std::min(dmin, PX2[i][j]); dmax = std::max(dmax, PX2[i][j]);
            cx += PX2[i][j]; cnt++;
        }
        cx /= cnt;
        double thick = BODY_W * BODY_H * f.girth * f.L * sc0;
        double len = dmax - dmin;
        if (len > 1e-6 && len < thick) {
            double k = thick / len;
            for (int i = 0; i <= NXG; ++i) for (int j = 0; j <= NYG; ++j)
                PX2[i][j] = cx + (PX2[i][j] - cx) * k;
        }
    }
    // 斜めを向くと頭と尾が画面で重なる。**奥の列から描かないと、
    // 奥にある尾が手前の頭を塗りつぶす**
    double dH = WXc[0] * WXc[0] + (WZc[0] + CAMD * NWATER) * (WZc[0] + CAMD * NWATER);
    double dT = WXc[NXG] * WXc[NXG] + (WZc[NXG] + CAMD * NWATER) * (WZc[NXG] + CAMD * NWATER);
    bool headFirst = (dH > dT);              // 頭のほうが奥なら頭から
    for (int ii = 0; ii < NXG; ++ii) {
        int i = headFirst ? ii : (NXG - 1 - ii);
        double u0 = (double)i / NXG, u1 = (double)(i + 1) / NXG;
        for (int j = 0; j < NYG; ++j) {
            double v0 = (double)j / NYG, v1 = (double)(j + 1) / NYG;
            double X4[4] = {PX2[i][j], PX2[i+1][j], PX2[i+1][j+1], PX2[i][j+1]};
            double Y4[4] = {PY2[i][j], PY2[i+1][j], PY2[i+1][j+1], PY2[i][j+1]};
            double U4[4] = {u0, u1, u1, u0}, V4[4] = {v0, v0, v1, v1};
            double S4[4] = {CSH[i] * RND[i][j],     CSH[i+1] * RND[i+1][j],
                            CSH[i+1] * RND[i+1][j+1], CSH[i] * RND[i][j+1]};
            double A4[4] = {ALP[i][j], ALP[i+1][j], ALP[i+1][j+1], ALP[i][j+1]};
            const int T1[3] = {0, 1, 2}, T2[3] = {0, 2, 3};
            for (int p = 0; p < 2; ++p) {
                const int* T = p ? T2 : T1;
                double xx[3], yy[3], uu[3], vv[3], ss[3], aa[3];
                for (int k = 0; k < 3; ++k) {
                    xx[k] = X4[T[k]]; yy[k] = Y4[T[k]]; uu[k] = U4[T[k]];
                    vv[k] = V4[T[k]]; ss[k] = S4[T[k]]; aa[k] = A4[T[k]];
                }
                warp_tri(xx, yy, uu, vv, ss, aa, fog);
            }
        }
    }
}
#endif

#if HAVE_PHOTO
static void draw_fish(const Fish& f) { draw_fish_full(f); }
#else
// 写真が無いときの手続き描画。形は出るが、板を並べた見え方になる
static void draw_fish(const Fish& f) {
    const int NS = 26;
    double cy = cos(f.yaw), sy = sin(f.yaw);
    double fx = cy * cos(f.pitch), fy = sin(f.pitch), fz = sy * cos(f.pitch);
    double ax = -sy, az = cy;

    double amp = 0.070 * f.L * f.ampS;   // 滑空中は 0 になり体がまっすぐになる
    double SX[NS + 1], SY[NS + 1], HH[NS + 1];
    for (int i = 0; i <= NS; ++i) {
        double s = (double)i / NS;
        double lat = body_wave(s * 0.9, f.phase, amp, f.bend);
        double along = f.L * (0.46 - s * 0.9);
        double sc;
        project(f.x + fx * along + ax * lat, f.y + fy * along, f.z + fz * along + az * lat,
                SX[i], SY[i], sc);
        HH[i] = 0.5 * body_height(s) * BODY_H * f.girth * f.L * sc;
    }

    // 正面を向いた魚が線に潰れないようにする。魚はリボン(板)ではなく厚みのある体なので、
    // 真正面から見ても「体の幅」ぶんは見えるはず。画面上の鼻先〜尾柄の長さが
    // 体の厚み(体高のおよそ0.4倍)を下回ったら、横方向だけそこまで広げる。
    {
        double cx = 0.0;
        for (int i = 0; i <= NS; ++i) cx += SX[i];
        cx /= (NS + 1);
        double dmin = SX[0], dmax = SX[0];
        for (int i = 1; i <= NS; ++i) { dmin = std::min(dmin, SX[i]); dmax = std::max(dmax, SX[i]); }
        double len = dmax - dmin;
        double thick = BODY_W * BODY_H * f.girth * f.L * (HH[6] > 0 ? (HH[6] / (0.5 * body_height(6.0 / NS) * BODY_H * f.girth * f.L)) : 1.0);
        if (len > 1e-6 && len < thick) {
            double k = thick / len;
            for (int i = 0; i <= NS; ++i) SX[i] = cx + (SX[i] - cx) * k;
        }
    }

    // 奥ほど水の色に沈む
    double dep = std::max(0.0, std::min(1.0, f.z / TD));
    double fog = 0.62 * dep * p_murk;
    // 構造色: 真横から見るほど強く光り、角度で青→緑に振れる
    double vx = f.x, vz = f.z + CAMD;
    double vl = sqrt(vx * vx + vz * vz); if (vl < 1e-6) vl = 1e-6;
    double side = fabs((ax * vx + az * vz) / vl);
    double irid = 0.42 + 0.58 * pow(side, 0.6);
    double hue  = 0.45 + 0.55 * side;

    // v(-1=背, +1=腹) と s から体色を作る。境目は smoothstep でぼかす
    // ネオンテトラの体側は上から順に「オリーブの背 / 青緑の構造色ライン /
    // 半透明の銀 / 後半だけ赤 / 白い腹」。境目はぼかしすぎると特徴が消えるので、
    // 幅を持たせつつコントラストは残す。
    // 胴体は CC0 写真から展開したテクスチャを (s, v) で引く。
    // 写真そのものなので鱗・虹色素胞・目・銀の腹がそのまま出る。
    // 角度による明るさ(真横から見るほど体側が見える)と水の減衰だけ後から掛ける。
    double shade = (0.78 + 0.30 * side) * f.tint;
    auto tex_color = [&](double v, double s, int& R, int& G, int& B) {
        double r, g, bl;
        // 写真が無いときの体色。金魚の朱は「背が濃く、腹に向かって薄くなって白に抜ける」。
        // v: -1=背 +1=腹
        double belly = smoothstep(0.05, 1.10, v);             // 腹側の白へなだらかに抜ける
        double dorsal = smoothstep(-1.00, -0.45, v);          // 背の濃い帯
        r  = 214 + 32 * dorsal + 34 * belly;
        g  =  58 + 40 * dorsal + 152 * belly;
        bl =  20 + 18 * dorsal + 140 * belly;
        // 更紗(赤白斑)。個体ごとに位置と量が変わる。境目は少しぼかす
        double n = sin(s * 6.1 + f.pat * 6.2831) + 0.7 * sin(s * 2.3 - v * 2.6 + f.pat * 3.7)
                 + 0.5 * sin(v * 3.4 + f.pat * 11.0);
        double white = smoothstep(0.15, 0.95, n * f.sara + (f.sara - 0.55));
        r  = r  * (1 - white) + 250 * white;
        g  = g  * (1 - white) + 238 * white;
        bl = bl * (1 - white) + 226 * white;
        // 上体のハイライト。丸い体に光が回っているのが分かる
        double hi = exp(-((v + 0.45) * (v + 0.45)) / 0.10) * 0.16;
        r += 255 * hi; g += 220 * hi; bl += 190 * hi;
        // 鱗。頭には無いので s で消す。細かすぎると縞に見えるので浅く
        double sc2 = 1.0 + 0.014 * sin(s * 118.0 + v * 3.0) * sin(v * 17.0)
                            * smoothstep(0.22, 0.38, s) * (1.0 - smoothstep(0.80, 0.98, s));
        r *= sc2; g *= sc2; bl *= sc2;
        // えらぶたの縁
        double gs = 0.205 + 0.030 * (1.0 - v * v);      // 中ほどが後ろへ張り出す弧
        double gill = smoothstep(gs - 0.022, gs, s) * (1.0 - smoothstep(gs, gs + 0.022, s));
        r *= 1.0 - 0.10 * gill; g *= 1.0 - 0.10 * gill; bl *= 1.0 - 0.10 * gill;
        // ---- 丸み。写真をそのまま貼ると板に見えるので、断面の丸さぶんの陰影を足す。
        // 体は横から見ると楕円なので、v(背-1〜腹+1) の位置で「こちらを向いている度合い」が
        // 決まる: nz = sqrt(1-v^2)。縁ほど光が回り込まず暗くなり、中ほどは明るい。
        // さらに体側の上寄りに、丸い面が光を返す帯(ハイライト)が出る
        double nz = sqrt(std::max(0.0, 1.0 - v * v));
        double round_ = 0.70 + 0.30 * nz;                    // 縁が落ちる
        double spec = exp(-((v + 0.30) * (v + 0.30)) / 0.055) * 0.20 * (0.5 + 0.5 * side);
        r = r * round_ + 255 * spec;
        g = g * round_ + 246 * spec;
        bl = bl * round_ + 232 * spec;
        r *= shade; g *= shade; bl *= shade;
        R = (int)(r * (1 - fog) + 14 * fog);
        G = (int)(g * (1 - fog) + 44 * fog);
        B = (int)(bl * (1 - fog) + 48 * fog);
    };
    (void)irid; (void)hue;

    // --- 胴体
    for (int i = 0; i < NS; ++i) {
        double x0 = SX[i], x1 = SX[i + 1];
        double y0 = SY[i], y1 = SY[i + 1];
        double h0 = HH[i], h1 = HH[i + 1];
        double s0 = (double)i / NS, s1 = (double)(i + 1) / NS;
        int xa = (int)floor(std::min(x0, x1)), xb = (int)ceil(std::max(x0, x1));
        if (xb < 0 || xa >= FW) continue;
        int n = std::max(1, xb - xa);
        for (int xi = xa; xi <= xb; ++xi) {
            if (xi < 0 || xi >= FW) continue;
            double u = (fabs(x1 - x0) < 0.5) ? 0.5 : (xi + 0.5 - x0) / (x1 - x0);
            if (u < 0) u = 0; if (u > 1) u = 1;
            double yc = y0 + (y1 - y0) * u, hh = h0 + (h1 - h0) * u, s = s0 + (s1 - s0) * u;
            double top = yc - hh, bot = yc + hh;
            int ya = (int)floor(top), yb = (int)ceil(bot);
            for (int yi = ya; yi <= yb; ++yi) {
                // 被覆率: この画素が体の内側に入っている割合
                double cov = std::min((double)yi + 1.0, bot) - std::max((double)yi, top);
                if (cov <= 0) continue;
                if (cov > 1) cov = 1;
                double v = (hh > 0.01) ? ((yi + 0.5) - yc) / hh : 0.0;
                if (v < -1) v = -1; if (v > 1) v = 1;
                int R, G, B; tex_color(v, s, R, G, B);
                blend_px(xi, yi, R, G, B, cov);
            }
        }
        (void)n;
    }

    // 画面に小さく写るときは、ひれを薄くする。実物のひれはほぼ透明で、
    // 遠い個体・引いた画角では体の輪郭しか見えない。鰭条は 0.75px の固定幅で
    // 描いているので、これをやらないと小さい個体ほど相対的に太くなって虫のように見える。
    double lenpx = sqrt((SX[NS] - SX[0]) * (SX[NS] - SX[0]) + (SY[NS] - SY[0]) * (SY[NS] - SY[0]));
    double memF = 0.10 + 0.90 * smoothstep(7.0, 30.0, lenpx);    // 膜
    double rayF = smoothstep(14.0, 40.0, lenpx);                 // 鰭条(細いので先に消す)

    // ひれは骨(鰭条)の間に薄い膜が張ったもの。ベタの三角形だと紙に見えるので、
    // 根元から先端へ向かう鰭条を細い三角形で何本か描き、その上に薄い膜を重ねる。
    // ひれの色は体色から取る。金魚のひれは体と同じ色素が薄く乗ったもので、
    // 朱の個体は朱いひれ、白勝ちの個体は白いひれになる(テトラのような青白ではない)
    int br_, bg_, bb_;
    tex_color(0.10, 0.86, br_, bg_, bb_);
    int fr  = (int)(br_ * 0.78 + 255 * 0.22);
    int fg2 = (int)(bg_ * 0.78 + 205 * 0.22);
    int fb  = (int)(bb_ * 0.78 + 165 * 0.22);
    // 根元(rx,ry)から out 方向へ len だけ伸び、spr 方向へ ±half に開く扇。
    // notch>0 で後縁の中央がえぐれる(尾びれの二叉)。鰭条を細い三角形で重ねる。
    auto fan = [&](double rx, double ry, double ox, double oy, double px2, double py2,
                   double len, double half, double notch, int nray, int memA, int rayA) {
        double ex[10], ey[10];
        int n = std::min(nray, 9);
        for (int k = 0; k <= n; ++k) {
            double q = (double)k / n * 2.0 - 1.0;
            double ext = len * ((1.0 - notch) + notch * fabs(q));
            ex[k] = rx + ox * ext + px2 * half * q;
            ey[k] = ry + oy * ext + py2 * half * q;
        }
        // 膜は根元ほど濃く、外縁で消える。実物のひれは縁に向かって薄くなるので、
        // 一様に塗ると紙を貼ったように見える。外側を薄く塗ってから内側を重ねる
        for (int k = 1; k <= n; ++k)
            tri(rx, ry, ex[k - 1], ey[k - 1], ex[k], ey[k], rgba(fr, fg2, fb, (int)(memA * 0.60 * memF)));
        for (int k = 1; k <= n; ++k) {
            const double IN = 0.55;
            tri(rx, ry, rx + (ex[k-1] - rx) * IN, ry + (ey[k-1] - ry) * IN,
                        rx + (ex[k]   - rx) * IN, ry + (ey[k]   - ry) * IN,
                rgba(fr, fg2, fb, (int)(memA * 0.75 * memF)));
        }
        // 鰭条は根元側だけ、うっすら。実物のひれはほぼ透明で、先端まで白い骨が通って
        // 見えることはない。全長を強い不透明度で描くと、ひれではなく針金の扇に見える
        // (とくに一番外の鰭条がひれの外縁と重なるので、縁が白い針金になる)
        if (rayA * rayF < 3.0) return;
        const double RAYLEN = 0.62;                          // 鰭条を描く長さ(ひれに対する割合)
        for (int k = 0; k <= n; ++k) {                       // 鰭条
            double vx2 = ex[k] - rx, vy2 = ey[k] - ry;
            double vl2 = sqrt(vx2 * vx2 + vy2 * vy2); if (vl2 < 1e-6) continue;
            double wx2 = -vy2 / vl2 * 0.75, wy2 = vx2 / vl2 * 0.75;
            tri(rx + wx2, ry + wy2, rx - wx2, ry - wy2,
                rx + vx2 * RAYLEN, ry + vy2 * RAYLEN, rgba(fr, fg2, fb, (int)(rayA * rayF)));
        }
    };

    // 尾のほうを向く画面上の単位ベクトル
    double uxs = SX[NS] - SX[NS - 4], uys = SY[NS] - SY[NS - 4];
    double ul2 = sqrt(uxs * uxs + uys * uys); if (ul2 < 1e-6) { uxs = 1; uys = 0; ul2 = 1; }
    uxs /= ul2; uys /= ul2;
    double sc0 = HH[6] / std::max(1e-9, 0.5 * body_height(6.0 / NS) * BODY_H * f.girth * f.L); // 画面倍率

    // --- 尾びれ。**体の波の続き**として、s=1 の先の中心線を出してから貼る。
    // 剛体の扇を尾柄に貼ると板を振っているようにしか見えない。金魚は尾が長いので、
    // 体が振り終わった後から尾が遅れて返ってくる。この遅れは進行波の位相差で自然に出る。
    {
        const int NF = 7;                        // 尾びれの中心線の分割
        const double SEND = 1.30;
        double FX[NF + 1], FY[NF + 1], NX[NF + 1], NY[NF + 1];
        for (int k = 0; k <= NF; ++k) {
            double sq = 0.930 + (SEND - 0.930) * (double)k / NF;   // 少し重ねて隙間を消す
            double q = std::max(0.0, (sq - 1.0) / std::max(1e-6, SEND - 1.0));   // 0=付け根 1=先
            // 付け根では体の波、そこから先は**ひれのしなり**へなめらかに移る
            double lat = body_wave(std::min(1.0, sq) * 0.9, f.phase, amp, f.bend) * (1.0 - q)
                       + fin_wave(q, f.phase, amp, 1.6, 1.8) * q;
            double along = f.L * (0.46 - sq * 0.9);
            double sc;
            project(f.x + fx * along + ax * lat, f.y + fy * along,
                    f.z + fz * along + az * lat, FX[k], FY[k], sc);
        }
        for (int k = 0; k <= NF; ++k) {
            int a2 = std::max(0, k - 1), b2 = std::min(NF, k + 1);
            double tx = FX[b2] - FX[a2], ty = FY[b2] - FY[a2];
            double tl = sqrt(tx * tx + ty * ty); if (tl < 1e-6) { tx = uxs; ty = uys; tl = 1; }
            NX[k] = -ty / tl; NY[k] = tx / tl;
        }
        // 上の葉と下の葉で遅れが違うので、尾は面としてねじれる。
        // 横から見たときに「ひらひら」して見えるのはこのねじれ
        double twist[NF + 1];
        for (int k = 0; k <= NF; ++k) {
            double q = (double)k / NF;
            twist[k] = 0.085 * f.L * sc0 * q * q * sin(f.phase * 0.85 - q * 2.2);
        }
        double UX[NF + 1], UY[NF + 1], LX[NF + 1], LY[NF + 1];
        for (int k = 0; k <= NF; ++k) {
            double q = (double)k / NF;
            double w = (0.042 + 0.150 * q * q) * f.L * sc0;
            UX[k] = FX[k] + NX[k] * w; UY[k] = FY[k] + NY[k] * w;
            LX[k] = FX[k] - NX[k] * w; LY[k] = FY[k] - NY[k] * w;
        }
        int memT = (int)(170 * memF), memT2 = (int)(110 * memF);
        for (int k = 0; k < NF; ++k) {
            uint32_t c = rgba(fr, fg2, fb, k < NF - 2 ? memT : memT2);
            tri(FX[k], FY[k], UX[k], UY[k], UX[k+1], UY[k+1], c);
            tri(FX[k], FY[k], UX[k+1], UY[k+1], FX[k+1], FY[k+1], c);
            tri(FX[k], FY[k], LX[k], LY[k], LX[k+1], LY[k+1], c);
            tri(FX[k], FY[k], LX[k+1], LY[k+1], FX[k+1], FY[k+1], c);
        }
    }

    // 写真が無いときは手続きで描く(形は出るが、板に見える)
    auto finstrip = [&](double s0, double s1, double sgn, double hMax, double lean,
                        int memA, int rayA) {
        const int NQ = 6;
        double RX[NQ + 1], RY[NQ + 1], TX[NQ + 1], TY[NQ + 1];
        for (int k = 0; k <= NQ; ++k) {
            double t = (double)k / NQ;
            double sq = s0 + (s1 - s0) * t;
            double fi = sq * NS; int i = clampi((int)fi, 0, NS - 1);
            double fr2 = fi - i;
            double rx = SX[i] + (SX[i + 1] - SX[i]) * fr2;
            double ry = SY[i] + (SY[i + 1] - SY[i]) * fr2;
            double hh = HH[i] + (HH[i + 1] - HH[i]) * fr2;
            RX[k] = rx; RY[k] = ry + sgn * hh * 0.94;
            double hgt = hMax * sin(M_PI * pow(t, 0.80));
            TX[k] = RX[k] + uxs * lean * hgt;
            TY[k] = RY[k] + sgn * hgt + uys * lean * hgt;
        }
        for (int k = 0; k < NQ; ++k) {
            uint32_t c = rgba(fr, fg2, fb, (int)(memA * memF));
            tri(RX[k], RY[k], RX[k+1], RY[k+1], TX[k+1], TY[k+1], c);
            tri(RX[k], RY[k], TX[k+1], TY[k+1], TX[k], TY[k], c);
        }
        if (rayA * rayF < 3.0) return;
        for (int k = 0; k <= NQ; ++k) {
            double vx2 = TX[k] - RX[k], vy2 = TY[k] - RY[k];
            double vl2 = sqrt(vx2 * vx2 + vy2 * vy2); if (vl2 < 1e-6) continue;
            double wx2 = -vy2 / vl2 * 0.6, wy2 = vx2 / vl2 * 0.6;
            tri(RX[k] + wx2, RY[k] + wy2, RX[k] - wx2, RY[k] - wy2,
                RX[k] + vx2 * 0.85, RY[k] + vy2 * 0.85, rgba(fr, fg2, fb, (int)(rayA * rayF)));
        }
    };
    finstrip(0.30, 0.62, -1.0, 0.175 * f.L * sc0, 0.55, 165, 26);
    finstrip(0.72, 0.86,  1.0, 0.105 * f.L * sc0, 0.60, 155, 24);
    {
        double sw = sin(f.pecPh * 0.7) * 0.10;
        finstrip(0.50, 0.60, 1.0, (0.085 + sw) * f.L * sc0, 0.75, 150, 22);
    }
    {
        int k0 = (int)(0.30 * NS);
        double sw = sin(f.pecPh);
        for (int sgn2 = 0; sgn2 < 2; ++sgn2) {
            double side2 = sgn2 ? 1.0 : 0.80;
            double ang = 0.55 + 0.42 * (sgn2 ? sw : -sw);
            double ox = uxs * ang, oy = 0.55 + uys * ang;
            double ol = sqrt(ox * ox + oy * oy); if (ol < 1e-6) ol = 1;
            fan(SX[k0], SY[k0] + HH[k0] * 0.34, ox / ol, oy / ol, uxs, uys,
                0.105 * f.L * sc0 * side2, 0.048 * f.L * sc0 * side2, 0.0, 5,
                (int)(150 * side2), (int)(24 * side2));
        }
    }

    // --- 目。**写真を貼っているときは描かない**(写真に目が写っているので二重になる)
    {
        double along = f.L * (0.46 - 0.085 * 0.9);
        double ex, ey, esc;
        // 金魚の目は頭の丸みの上のほうに付く(テトラより高く、少し大きい)
        project(f.x + fx * along, f.y + fy * along + 0.012 * f.L, f.z + fz * along, ex, ey, esc);
        double er = std::max(1.0, 0.040 * f.L * esc);
        // 虹彩(銀に光る)→瞳(黒)→ハイライト
        olivec_circle(OC, (int)ex, (int)ey, (int)(er + 0.5),
                      rgb((int)((186 + 60 * side) * (1 - fog) + 14 * fog),
                          (int)((132 + 55 * side) * (1 - fog) + 44 * fog),
                          (int)(( 52 + 40 * side) * (1 - fog) + 48 * fog)));
        olivec_circle(OC, (int)ex, (int)ey, std::max(1, (int)(er * 0.66)),
                      rgb((int)(16 * (1 - fog) + 14 * fog), (int)(18 * (1 - fog) + 44 * fog),
                          (int)(22 * (1 - fog) + 48 * fog)));
        if (er > 2.0)
            olivec_circle(OC, (int)(ex - er * 0.28), (int)(ey - er * 0.30), std::max(1, (int)(er * 0.3)),
                          rgb(225, 233, 238));
    }
}

#endif

// ---------------------------------------------------------------- 近傍探索
// 群れの計算は「見える距離(Rsee)にいる仲間」だけを見る。総当たりだと N² で、
// 900匹だと 81万組/サブステップになって間に合わない(実測 400匹で 4.9ms/フレーム)。
// 一辺 Rsee の格子に入れておけば、自分と隣の 27 セルを見るだけで漏れなく拾える。
static const double RSEE = 0.26;      // ここまでが仲間として見える(体が大きいので広い)
static const double RSEP = 0.135;     // これより近いと離れる(体長の1.3倍ほど)
static int GX = 1, GY = 1, GZ = 1;
static std::vector<int> gStart, gFill, gItem;

static inline int grid_cell(const Fish& f) {
    int ix = clampi((int)((f.x + TW * 0.5) / RSEE), 0, GX - 1);
    int iy = clampi((int)(f.y / RSEE), 0, GY - 1);
    int iz = clampi((int)(f.z / RSEE), 0, GZ - 1);
    return (iz * GY + iy) * GX + ix;
}
static void grid_init() {
    GX = std::max(1, (int)(TW / RSEE) + 1);
    GY = std::max(1, (int)(TH / RSEE) + 1);
    GZ = std::max(1, (int)(TD / RSEE) + 1);
}
static void grid_build() {
    int nc = GX * GY * GZ;
    gStart.assign(nc + 1, 0);
    gFill.assign(nc, 0);
    gItem.resize(fish.size());
    for (auto& f : fish) gStart[grid_cell(f) + 1]++;
    for (int c = 0; c < nc; ++c) gStart[c + 1] += gStart[c];
    for (size_t i = 0; i < fish.size(); ++i) {
        int c = grid_cell(fish[i]);
        gItem[gStart[c] + gFill[c]++] = (int)i;
    }
}

// ---------------------------------------------------------------- 生成
static Fish make_fish() {
        Fish f;
        f.x = rnds() * TW * 0.46;
        f.y = TH * (0.06 + 0.88 * rnd());
        f.z = TD * (0.12 + 0.76 * rnd());
        f.yaw = rnd() * 6.2831853;
        f.pitch = rnds() * 0.15;
        // 和金は 8〜13cm(尾を除く胴の長さ)。当歳魚の小さいのが混じる
        double sz = (rnd() < 0.22) ? (0.55 + 0.20 * rnd())     // 当歳魚
                                   : (0.88 + 0.26 * rnd());    // 親魚
        f.L = 0.105 * sz;
        f.girth = 0.88 + 0.24 * rnd();
        if (sz < 0.80) f.girth *= 0.92;                        // 小さい個体はほっそり
        f.speed = 0.03;
        f.phase = rnd() * 6.2831853;
        f.freq = 1.6;
        f.bend = 0.0;
        f.dartT = 0.0;
        f.gait = (rnd() < 0.5) ? 0 : 1;
        f.beats = 2.0 + 2.0 * rnd();
        f.coastT = 0.3 + 0.7 * rnd();
        f.ampS = 1.0;
        f.cruiseK = 0.78 + 0.50 * rnd();
        f.tint = (float)(0.9 + 0.2 * rnd());
        // 模様の個体差。素赤(全身朱)から白勝ちの更紗まで
        f.pat  = (float)rnd();
        double q = rnd();
        f.sara = (float)(q < 0.35 ? 0.05 + 0.12 * rnd()          // 素赤(全身朱)
                       : q < 0.85 ? 0.30 + 0.30 * rnd()          // 赤勝ちの更紗
                                  : 0.70 + 0.25 * rnd());        // 白勝ち(少なめ)
        f.pecPh = rnd() * 6.2831853;
        return f;
}

// 密度から匹数を出す
static int target_count() {
    double liters = TW * TH * TD * 1000.0;
    return clampi((int)(p_dens * liters + 0.5), 1, (int)MAX_FISH);
}

// 水槽が変わった/密度が変わったときに匹数を合わせる。
// 縮んだら外に出た個体を落とし、広がったら足りない分を足す。
// 全部作り直すと、スライダを動かすたびに群れが消えて作り直しになってしまう
static void refit_fish() {
    for (size_t i = 0; i < fish.size(); ) {
        const Fish& f = fish[i];
        if (fabs(f.x) > TW * 0.5 || f.y < 0.0 || f.y > TH) { fish[i] = fish.back(); fish.pop_back(); }
        else ++i;
    }
    int want = target_count();
    while ((int)fish.size() > want) fish.pop_back();
    while ((int)fish.size() < want) fish.push_back(make_fish());
}

static void spawn_fish() {
    fish.clear();
    int n = target_count();
    fish.reserve(n);
    for (int i = 0; i < n; ++i) fish.push_back(make_fish());
}

// ---------------------------------------------------------------- ABI
extern "C" {

KEEP int sim_w() { return FW; }
KEEP int sim_h() { return FH; }

// 画角かカメラ距離が変わったら、水槽の大きさもそれに合わせて作り直す
static void apply_view() {
    FOC = FW * CAMD / VIEWW;
    tank_from_view();
    grid_init();
    build_static_bg();
    refit_fish();
}

KEEP void sim_reset() {
    FOC = FW * CAMD / VIEWW;
    tank_from_view();
    grid_init();
    build_static_bg();
    spawn_fish();
    foods.clear();
    startle = 0.0; simTime = 0; frameNo = 0;
}

KEEP int sim_init(int, int) {
    px.assign((size_t)FW * FH, 0);
    build_upsample_tables();
    sim_reset();
    return 1;
}

KEEP void sim_set(int id, double v) {
    switch (id) {
    // 0 は密度[匹/L]。水槽の大きさが画角で変わるので、匹数ではなく密度で持つ
    case 0: if (fabs(v - p_dens) > 1e-9) { p_dens = v; refit_fish(); } break;
    case 1: p_school = v; break;
    case 2: p_speed  = v; break;
    case 3: if (fabs(v - p_murk) > 1e-9) { p_murk = v; build_static_bg(); } break;
    case 4: p_light  = v; break;
    case 6: p_hud    = v; break;
    case 20: {                                  // モデル確認モード
        bool on = v > 0.5;
        if (on != p_model) {
            p_model = on;
            if (on) {
                fish.clear();
                Fish f = make_fish();
                f.L = 0.105; f.girth = 1.0; f.tint = 1.0f;
                f.sara = 0.55f; f.pat = 0.30f;
                f.yaw = M_PI; f.pitch = 0.0;     // 左を向く
                f.speed = f.L * 0.8;
                fish.push_back(f);
            } else refit_fish();
        }
        break; }
    // 画角とカメラ距離を変えると、水槽の大きさも作り直す(隅が見えない大きさ + WING)。
    // 範囲を切るのは水槽の都合ではなく、匹数の上限(MAX_FISH)で密度が保てなくなるから
    case 7: { double w = std::max(0.10, std::min(VIEWW_MAX, v));
              if (fabs(w - VIEWW) > 1e-9) { VIEWW = w; apply_view(); } break; }
    case 8: { double d = std::max(CAMD_MIN, std::min(5.0, v));
              if (fabs(d - CAMD) > 1e-9) { CAMD = d; apply_view(); } break; }
    }
}

KEEP double sim_get(int id) {
    switch (id) {
    case 0: return (double)fish.size();
    case 2: return (double)foods.size();
    case 6: { int n = 0; for (auto& f : fish) if (f.gait == 1) n++;
              return fish.empty() ? 0 : (double)n / fish.size(); }   // 滑空中の割合
    case 3: return simTime;
    case 4: { double a = 0; for (auto& f : fish) a += f.freq * f.ampS; return fish.empty() ? 0 : a / fish.size(); }
    case 5: { double a = 0; for (auto& f : fish) a += f.speed / f.L; return fish.empty() ? 0 : a / fish.size(); }
    case 7: return TW * TH * TD * 1000.0;    // いまの水槽の容量 [L] (画角で変わる)
    case 8: return TW;                       // 幅 [m]
    case 9: return TH;                       // 高さ [m]
    }
    return 0.0;
}

KEEP void sim_action(int id) {
    // 餌は「見えているところの上のほう」に落とす。水面は画面の外(背の高い水槽なので)、
    // そこから落とすと 0.02m/s では画面に入るまで何十秒もかかってしまう。
    // 画面のちょうど上端でもなく少し内側にするのは、集まってくるところが画面に入るように
    if (id == 0) { Food fd;
                   fd.z = TD * (0.3 + 0.4 * rnd());
                   fd.x = rnds() * std::min(TW * 0.45, view_half_w(fd.z) * 0.85);
                   fd.y = std::min(TH * 0.97, CAMY + view_half_h(fd.z) * 0.55);
                   fd.t = 0;
                   if (foods.size() < 40) foods.push_back(fd); }
    else if (id == 1) sim_reset();
    else if (id == 2) startle = 1.0;
}

// クリックしたところに餌を落とす。水面は画面の外なので、そこから落としても届かない。
// 投影と同じ屈折(z/n + CAMD)で逆算する
KEEP void sim_click(double nx, double ny) {
    double zc = TD * 0.5, sc = FOC / (zc / NWATER + CAMD);
    Food fd;
    fd.x = (nx * FW - FW * 0.5) / sc;
    fd.y = CAMY - (ny * FH - FH * 0.5) / sc;
    fd.z = zc + rnds() * TD * 0.2;
    fd.t = 0;
    fd.x = std::max(-TW * 0.45, std::min(TW * 0.45, fd.x));
    fd.y = std::max(0.02, std::min(TH * 0.97, fd.y));
    if (foods.size() < 40) foods.push_back(fd);
}

KEEP void sim_step(int frames) {
    for (int fr = 0; fr < frames; ++fr) {
        for (int ss = 0; ss < SUBSTEPS; ++ss) {
            const double dt = DT;
            simTime += dt;
            if (startle > 0) startle = std::max(0.0, startle - dt * 0.7);

            for (size_t i = 0; i < foods.size();) {
                foods[i].y -= 0.02 * dt; foods[i].t += dt;
                if (foods[i].y < 0.012 || foods[i].t > 25.0) { foods[i] = foods.back(); foods.pop_back(); }
                else ++i;
            }

            size_t N = fish.size();
            grid_build();
            if (p_model && !fish.empty()) {
                // 画面中央に置き直す。前へ進む力はそのままなので、体もひれも普通に動く
                Fish& f = fish[0];
                f.x = 0.0; f.y = CAMY; f.z = TD * 0.5;
                f.yaw = M_PI; f.pitch = 0.0; f.bend = 0.0;
                f.gait = 0; f.ampS = 1.0;
                f.speed = f.L * 0.8 * p_speed;
                f.phase += 2.0 * M_PI * f.freq * dt;
                f.pecPh += 3.0 * dt;
                f.freq = std::max(0.8, std::min(6.0, 0.30 * f.speed / (2.0 * 0.070 * f.L)));
                continue;
            }
            for (size_t i = 0; i < N; ++i) {
                Fish& f = fish[i];
                double sepx = 0, sepy = 0, sepz = 0, alx = 0, aly = 0, alz = 0, cox = 0, coy = 0, coz = 0;
                int nal = 0, nco = 0;
                const double Rsep = RSEP, Rsee = RSEE;
                int ix = clampi((int)((f.x + TW * 0.5) / RSEE), 0, GX - 1);
                int iy = clampi((int)(f.y / RSEE), 0, GY - 1);
                int iz = clampi((int)(f.z / RSEE), 0, GZ - 1);
                for (int jz = std::max(0, iz - 1); jz <= std::min(GZ - 1, iz + 1); ++jz)
                for (int jy = std::max(0, iy - 1); jy <= std::min(GY - 1, iy + 1); ++jy)
                for (int jx = std::max(0, ix - 1); jx <= std::min(GX - 1, ix + 1); ++jx) {
                    int c = (jz * GY + jy) * GX + jx;
                    for (int q = gStart[c]; q < gStart[c + 1]; ++q) {
                    size_t k = (size_t)gItem[q];
                    if (k == i) continue;
                    double dx = fish[k].x - f.x, dy = fish[k].y - f.y, dz = fish[k].z - f.z;
                    double d2 = dx * dx + dy * dy + dz * dz;
                    if (d2 > Rsee * Rsee || d2 < 1e-9) continue;
                    double d = sqrt(d2);
                    if (d < Rsep) { double w = (Rsep - d) / Rsep / d; sepx -= dx * w; sepy -= dy * w; sepz -= dz * w; }
                    alx += cos(fish[k].yaw); aly += sin(fish[k].pitch); alz += sin(fish[k].yaw); nal++;
                    cox += dx; coy += dy; coz += dz; nco++;
                    }
                }
                double dx = 0, dy = 0, dz = 0, sw = p_school;
                dx += sepx * 5.5; dy += sepy * 5.5; dz += sepz * 5.5;
                // 金魚は群れを作らない。互いに避けはするが、向きも位置も揃えようとしない
                // (テトラは 1.4 / 0.9。同じ値にすると金魚が群泳して別の魚に見える)
                if (nal) { dx += alx / nal * 0.30 * sw; dy += aly / nal * 0.30 * sw; dz += alz / nal * 0.30 * sw; }
                if (nco) { dx += cox / nco * 0.18 * sw; dy += coy / nco * 0.18 * sw; dz += coz / nco * 0.18 * sw; }

                const double m = 0.05;
                double lx = -TW * 0.5 + m, rx = TW * 0.5 - m;
                if (f.x < lx) dx += (lx - f.x) * 70.0;
                if (f.x > rx) dx -= (f.x - rx) * 70.0;
                if (f.y < m)      dy += (m - f.y) * 80.0;
                if (f.y > TH - m) dy -= (f.y - (TH - m)) * 80.0;
                if (f.z < m * 0.7)      dz += (m * 0.7 - f.z) * 70.0;
                if (f.z > TD - m * 0.7) dz -= (f.z - (TD - m * 0.7)) * 70.0;

                double hunger = 0.0;
                for (size_t q = 0; q < foods.size();) {
                    Food& fd = foods[q];
                    double ex = fd.x - f.x, ey = fd.y - f.y, ez = fd.z - f.z;
                    double e2 = ex * ex + ey * ey + ez * ez;
                    // 口が届いたら食べる。餌はそこで消える(底まで落ちない)
                    double eat = f.L * 0.28;
                    if (e2 < eat * eat) {
                        foods[q] = foods.back(); foods.pop_back();
                        f.dartT = 0.0;
                        continue;
                    }
                    if (e2 < 0.09) {
                        double e = sqrt(e2) + 1e-6, w = 5.0 / (e + 0.05);
                        dx += ex / e * w; dy += ey / e * w; dz += ez / e * w;
                        hunger = std::max(hunger, 1.0 - e / 0.3);
                    }
                    ++q;
                }
                if (startle > 0.01) {
                    dx += f.x * 12.0 * startle; dz += (f.z - TD * 0.5) * 12.0 * startle;
                    dy += 4.0 * startle * rnds();
                }
                dx += rnds() * 1.3; dy += rnds() * 0.7; dz += rnds() * 1.3;

                double dl = sqrt(dx * dx + dy * dy + dz * dz);
                if (dl > 1e-6) {
                    double tYaw = atan2(dz, dx);
                    double tPit = asin(std::max(-1.0, std::min(1.0, dy / dl)));
                    double dyaw = tYaw - f.yaw;
                    while (dyaw >  M_PI) dyaw -= 2 * M_PI;
                    while (dyaw < -M_PI) dyaw += 2 * M_PI;
                    // 金魚は体が大きくて重いので、テトラのようには曲がれない。
                    // そのぶん曲がるときは体を大きく弓なりにする(見ていて金魚らしい所)
                    const double maxRate = 1.7;
                    double rate = std::max(-maxRate, std::min(maxRate, dyaw * 1.8));
                    f.yaw += rate * dt;
                    double dpit = std::max(-0.6, std::min(0.6, tPit * 0.6 - f.pitch));
                    f.pitch += dpit * 1.5 * dt;
                    f.pitch = std::max(-0.36, std::min(0.36, f.pitch));
                    double bendT = -rate / maxRate * 0.55 * f.L;   // 金魚は曲がるとき体を大きく弓なりにする
                    f.bend += (bendT - f.bend) * 5.0 * dt;
                }

                if (hunger > 0.5 && rnd() < 0.02) f.dartT = 0.35;
                if (f.dartT > 0) f.dartT -= dt;
                // 巡航はおよそ 0.8 体長/秒。金魚は水槽の中ではゆったり漂うように泳ぐ
                // (テトラは 1.5 体長/秒。同じ値にすると金魚が忙しなく見える)
                double cruise = f.L * 0.8 * f.cruiseK * p_speed;
                double urge = 1.0 + 1.1 * hunger + 2.4 * startle + (f.dartT > 0 ? 2.0 : 0.0);

                // --- burst-and-coast
                // 尾を数回振って加速し、体をまっすぐにして惰性で滑る。滑空中は二次抗力で減速。
                // 金魚はテトラほど極端に滑らない(体が大きく、抗力係数 Cd は体長に反比例するので
                // 惰性がよく効く)。振る回数を増やし、滑空を短めにしてある。
                const double GLIDE_CD = 3.2;               // [1/m] 滑空時の抗力係数
                if (f.gait == 0) {
                    f.ampS += (1.0 - f.ampS) * 14.0 * dt;
                    f.speed += (cruise * urge * 1.30 - f.speed) * 5.0 * dt;
                    f.beats -= f.freq * dt;                // 振った回数を数える
                    if (f.beats <= 0.0) {
                        f.gait = 1;
                        f.coastT = (0.35 + 1.10 * rnd()) / std::max(0.35, urge);
                    }
                } else {
                    f.ampS += (0.0 - f.ampS) * 11.0 * dt;  // 体をまっすぐに戻す
                    f.speed -= GLIDE_CD * f.speed * f.speed * dt;
                    f.coastT -= dt;
                    // 遅くなりすぎたか、滑空の時間が尽きたら、また数回振る
                    if (f.coastT <= 0.0 || f.speed < cruise * urge * 0.55) {
                        f.gait = 0;
                        f.beats = 3.0 + 3.0 * rnd();
                    }
                }
                if (f.speed < 0.004) f.speed = 0.004;

                double cy = cos(f.yaw), sy = sin(f.yaw), cp = cos(f.pitch);
                f.x += cy * cp * f.speed * dt;
                f.y += sin(f.pitch) * f.speed * dt;
                f.z += sy * cp * f.speed * dt;
                f.x = std::max(-TW * 0.5 + 0.006, std::min(TW * 0.5 - 0.006, f.x));
                f.y = std::max(0.013, std::min(TH - 0.009, f.y));
                f.z = std::max(0.013, std::min(TD - 0.009, f.z));

                // 胸びれは止まっているときほど忙しなく漕ぐ(ホバリング)。
                // 泳いでいるときは畳み気味なので、ゆっくりになる
                f.pecPh += (2.2 + 3.4 * std::max(0.0, 1.0 - f.speed / (cruise + 1e-9))) * dt;
                if (f.pecPh > 1e6) f.pecPh = fmod(f.pecPh, 6.2831853);

                // 尾の振り数はストローハル数から: f = St U / (2A)
                double A = 0.070 * f.L;
                f.freq = std::max(0.8, std::min(14.0, 0.30 * f.speed / (2.0 * A)));
                f.phase += 2.0 * M_PI * f.freq * dt;
                if (f.phase > 1e6) f.phase = fmod(f.phase, 2.0 * M_PI);
            }

        }
        frameNo++;
    }
}

KEEP uint8_t* sim_render() {
    OC = olivec_canvas(px.data(), FW, FH, FW);
    draw_water();

    // --- 魚は奥から手前へ
    static std::vector<int> order;
    order.resize(fish.size());
    for (size_t i = 0; i < fish.size(); ++i) order[i] = (int)i;
    std::sort(order.begin(), order.end(), [](int a, int b) { return fish[a].z > fish[b].z; });
    for (int idx : order) draw_fish(fish[idx]);

    for (auto& fd : foods) {
        double sxp, syp, scp;
        project(fd.x, fd.y, fd.z, sxp, syp, scp);
        olivec_circle(OC, (int)sxp, (int)syp, std::max(1, (int)(0.0012 * scp)), rgb(196, 148, 72));
    }
    if (p_hud >= 0.5) {
        char buf[160];
        snprintf(buf, sizeof(buf), "WAKIN  n=%d   tail %.1f Hz   speed %.1f BL/s",
                 (int)fish.size(), sim_get(4), sim_get(5));
        olivec_text(OC, buf, 14, 12, olivec_default_font, 2, rgb(150, 220, 225));
    }
    return (uint8_t*)px.data();
}

}  // extern "C"

// ---------------------------------------------------------------- native self-test
#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <cstdio>
#include <cstdlib>
int main(int argc, char** argv) {
    sim_init(0, 0);
    int steps = argc > 1 ? atoi(argv[1]) : 400;
    const char* out = argc > 2 ? argv[2] : "suisou_preview.png";
    if (argc > 3) sim_set(0, atof(argv[3]));     // 4番目 = 密度[匹/L] (匹数は水槽の大きさで決まる)
    if (argc > 4) sim_set(6, atof(argv[4]));     // 5番目 = 0 で HUD なし
    int feedAt = argc > 5 ? atoi(argv[5]) : 0;   // 6番目 = このフレームで餌を撒く(0でなし)
    if (argc > 6) sim_set(7, atof(argv[6]));     // 7番目 = 画角[m]
    if (argc > 8 && atof(argv[8]) > 0.5) sim_set(20, 1);   // 9番目 = 1 でモデル確認モード
    if (argc > 7) sim_set(8, atof(argv[7]));     // 8番目 = 遠近感(カメラ距離[m])
    for (int i = 0; i < steps; ++i) {
        if (feedAt && i == feedAt) for (int k = 0; k < 6; ++k) sim_action(0);
        sim_step(1); sim_render();
    }
    uint8_t* p = sim_render();
    printf("suisou_os native: %dx%d frames=%d fish=%d\n",
           FW, FH, steps, (int)fish.size());
    printf("  tail %.2f Hz   speed %.2f BL/s   (どちらも積分の結果)\n", sim_get(4), sim_get(5));
    stbi_write_png(out, FW, FH, 4, p, FW * 4);
    return 0;
}
#endif
