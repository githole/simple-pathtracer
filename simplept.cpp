// Copyright (c) 2014 hole
// This software is released under the MIT License (http://kagamin.net/hole/license.txt).
// A part of this software is based on smallpt (http://www.kevinbeason.com/smallpt/) and
// released under the MIT License (http://kagamin.net/hole/smallpt-license.txt).
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>

const double PI = 3.14159265358979323846;
const double INF = 1e20;
const double EPS = 1e-5;
const double MaxDepth = 5;

// *** その他の関数 ***
inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; } 
inline int toInt(double x){ return int(pow(clamp(x),1/2.2)*255+.5); } 
inline double rand01() { return (double)rand()/RAND_MAX; }

// *** データ構造 ***
struct Vec {
	double x, y, z;
	Vec(const double x_ = 0, const double y_ = 0, const double z_ = 0) : x(x_), y(y_), z(z_) {}
	inline Vec operator+(const Vec &b) const {return Vec(x + b.x, y + b.y, z + b.z);}
	inline Vec operator-(const Vec &b) const {return Vec(x - b.x, y - b.y, z - b.z);}
	inline Vec operator*(const double b) const {return Vec(x * b, y * b, z * b);}
	inline Vec operator/(const double b) const {return Vec(x / b, y / b, z / b);}
	inline const double LengthSquared() const { return x*x + y*y + z*z; }
	inline const double Length() const { return sqrt(LengthSquared()); }
};
inline Vec operator*(double f, const Vec &v) { return v * f; }
inline Vec Normalize(const Vec &v) { return v / v.Length(); }
// 要素ごとの積をとる
inline const Vec Multiply(const Vec &v1, const Vec &v2) {
	return Vec(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
}
inline const double Dot(const Vec &v1, const Vec &v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}
inline const Vec Cross(const Vec &v1, const Vec &v2) {
	return Vec((v1.y * v2.z) - (v1.z * v2.y), (v1.z * v2.x) - (v1.x * v2.z), (v1.x * v2.y) - (v1.y * v2.x));
}

// RGB
typedef Vec Color;
static const Color BackgroundColor(0.0, 0.0, 0.0);

// スペクトル
typedef double Spectrum;
static const Spectrum BackgroundSpectrum = 0.0;

struct Ray {
	Vec org, dir;
	Ray(const Vec org_, const Vec &dir_) : org(org_), dir(dir_) {}
};

enum ReflectionType {
	DIFFUSE,    // 完全拡散面。いわゆるLambertian面。
	SPECULAR,   // 理想的な鏡面。
	REFRACTION, // 理想的なガラス的物質。
};

struct Sphere {
	double radius;
	Vec position;
	Spectrum emission, color;
	ReflectionType ref_type;

	Sphere(const double radius_, const Vec &position_, const Spectrum &emission_, const Spectrum &color_, const ReflectionType ref_type_) :
	  radius(radius_), position(position_), emission(emission_), color(color_), ref_type(ref_type_) {}
	// 入力のrayに対する交差点までの距離を返す。交差しなかったら0を返す。
	const double intersect(const Ray &ray) {
		Vec o_p = position - ray.org;
		const double b = Dot(o_p, ray.dir), det = b * b - Dot(o_p, o_p) + radius * radius;
		if (det >= 0.0) {
			const double sqrt_det = sqrt(det);
			const double t1 = b - sqrt_det, t2 = b + sqrt_det;
			if (t1 > EPS)		return t1;
			else if(t2 > EPS)	return t2;
		}
		return 0.0;
	}
};

Sphere spheres[] = {
	Sphere(1.0, Vec(36, 75, 45.0),     Spectrum(),  Spectrum(), DIFFUSE),//照明
	Sphere(1e5, Vec( 1e5+1,40.8,81.6), Spectrum(),  Spectrum(0.75),DIFFUSE),// 左
	Sphere(1e5, Vec(-1e5+99,40.8,81.6),Spectrum(),  Spectrum(0.75),DIFFUSE),// 右Z
	Sphere(1e5, Vec(50,40.8, 1e5),     Spectrum(),  Spectrum(0.75),DIFFUSE),// 奥
	Sphere(1e5, Vec(50,40.8,-1e5+170), Spectrum(),  Spectrum(), DIFFUSE),// 手前
	Sphere(1e5, Vec(50, 1e5, 81.6),    Spectrum(),  Spectrum(0.75),DIFFUSE),// 床
	Sphere(1e5, Vec(50,-1e5+81.6,81.6),Spectrum(),  Spectrum(0.75),DIFFUSE),// 天井
	Sphere(30, Vec(50, 40, 58),        Spectrum(),  Spectrum(0.99), REFRACTION),//ガラス
};
const int LightID = 0;

inline Spectrum light_intensity(const double wavelength) {
	return (120.0 / 95.0);
//	return (12.0/95.0) * (exp(-pow(wavelength - 600, 2)/10000) * 450 + 150) / 650.0;
}


// *** レンダリング用関数 ***
// シーンとの交差判定関数
inline bool intersect_scene(const Ray &ray, double *t, int *id) {
	const double n = sizeof(spheres) / sizeof(Sphere);
	*t  = INF;
	*id = -1;
	for (int i = 0; i < int(n); i ++) {
		double d = spheres[i].intersect(ray);
		if (d > 0.0 && d < *t) {
			*t  = d;
			*id = i;
		}
	}
	return *t < INF;
}


// 光源上の点をサンプリングして直接光を計算する
Spectrum direct_radiance_sample(const Vec &v0, const Vec &normal, const int id,  const double wavelength) {
	// 光源上の一点をサンプリングする
	const double r1 = 2 * PI * rand01();
	const double r2 = 1.0 - 2.0 * rand01();
	const Vec light_pos = spheres[LightID].position + ((spheres[LightID].radius + EPS) * Vec(sqrt(1.0 - r2*r2) * cos(r1), sqrt(1.0 - r2*r2) * sin(r1), r2));
	
	// サンプリングした点から計算
	const Vec light_normal = Normalize(light_pos - spheres[LightID].position);
	const Vec light_dir = Normalize(light_pos - v0);
	const double dist2 = (light_pos - v0).LengthSquared();
	const double dot0 = Dot(normal, light_dir);
	const double dot1 = Dot(light_normal, -1.0 * light_dir);

	if (dot0 >= 0 && dot1 >= 0) {
		const double G = dot0 * dot1 / dist2;
		double t; // レイからシーンの交差 位置までの距離
		int id_; // 交差したシーン内オブジェクトのID
		intersect_scene(Ray(v0, light_dir), &t, &id_);
		if (fabs(sqrt(dist2) - t) < 1e-3) {	
			const Spectrum emission = light_intensity(wavelength);
			return emission * spheres[id].color * (1.0 / PI) * G / (1.0 / (4.0 * PI * pow(spheres[LightID].radius, 2.0)));
		}
	}
	return Spectrum(0.0);
}

// ray方向からの放射輝度を求める
Spectrum radiance(const Ray &ray, const int depth, const double wavelength) {
	double t; // レイからシーンの交差位置までの距離
	int id;   // 交差したシーン内オブジェクトのID
	if (!intersect_scene(ray, &t, &id)) {
		return BackgroundSpectrum;
	}
	const Sphere &obj = spheres[id];
	const Vec hitpoint = ray.org + t * ray.dir; // 交差位置
	const Vec normal  =  Normalize(hitpoint - obj.position); // 交差位置の法線
	Vec orienting_normal = Dot(normal, ray.dir) < 0.0 ? normal : (-1.0 * normal); // 交差位置の法線（物体からのレイの入出を考慮）
	// 色の反射率最大のものを得る。ロシアンルーレットで使う。
	// ロシアンルーレットの閾値は任意だが色の反射率等を使うとより良い。
	double russian_roulette_probability = obj.color;
	double hoge = rand01();
	// 一定以上レイを追跡したらロシアンルーレットを実行し追跡を打ち切るかどうかを判断する
	if (depth > MaxDepth) {
		if (hoge >= russian_roulette_probability) {
			return Spectrum(0.0);
		}
	} else
		russian_roulette_probability = 1.0; // ロシアンルーレット実行しなかった
	

	switch (obj.ref_type) {
	case DIFFUSE: {
		if (id != LightID) {
			const int shadow_ray = 1;
			Spectrum direct_light(0.0);
			for (int i = 0; i < shadow_ray; i ++) {
				direct_light = direct_light + direct_radiance_sample(hitpoint, orienting_normal, id, wavelength) / shadow_ray;
			}

			// orienting_normalの方向を基準とした正規直交基底(w, u, v)を作る。この基底に対する半球内で次のレイを飛ばす。
			Vec w, u, v;
			w = orienting_normal;
			if (fabs(w.x) > 0.1)
				u = Normalize(Cross(Vec(0.0, 1.0, 0.0), w));
			else
				u = Normalize(Cross(Vec(1.0, 0.0, 0.0), w));
			v = Cross(w, u);
			// コサイン項を使った重点的サンプリング
			const double r1 = 2 * PI * rand01();
			const double r2 = rand01(), r2s = sqrt(r2);
			Vec dir = Normalize((u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1.0 - r2)));
			return (direct_light + obj.color * radiance(Ray(hitpoint, dir), depth+1, wavelength)) / russian_roulette_probability;
		} else if (depth == 0) {
			return light_intensity(wavelength);
		} else
			return Spectrum();

	} break;
	case SPECULAR: {
		// 完全鏡面なのでレイの反射方向は決定的。
		// ロシアンルーレットの確率で除算するのは上と同じ。
		// 直接光サンプリングする
		double lt;
		int lid;
		Ray reflection_ray = Ray(hitpoint, ray.dir - normal * 2.0 * Dot(normal, ray.dir));
		intersect_scene(reflection_ray, &lt, &lid);
		Spectrum direct_light(0.0);
		if (lid == LightID)
			direct_light = light_intensity(wavelength);

		return (direct_light + obj.color * radiance(reflection_ray, depth+1, wavelength)) / russian_roulette_probability;
	} break;
	case REFRACTION: {
		Ray reflection_ray = Ray(hitpoint, ray.dir - normal * 2.0 * Dot(normal, ray.dir));

		// 反射方向からの直接光サンプリングする
		double lt;
		int lid;
		intersect_scene(reflection_ray, &lt, &lid);
		Spectrum direct_light(0.0);
		if (lid == LightID)
			direct_light = light_intensity(wavelength);

		bool into = Dot(normal, orienting_normal) > 0.0; // レイがオブジェクトから出るのか、入るのか

		// Snellの法則
		const double nc = 1.0; // 真空の屈折率

		// http://refractiveindex.info/?group=CDGM&material=H-ZF62
		// CDGM社のH-ZF62という光学ガラス。屈折率の変化が大きい
		const double nt = 
			sqrt(3.5139564 - 2.4812508E-2 * pow(wavelength, 2) +
				4.6252559E-2 * pow(wavelength, -2) +
				9.1313596E-3 * pow(wavelength, -4) -
				1.0777108E-3 * pow(wavelength, -6) +
				1.0819677E-4 * pow(wavelength, -8));

		// ダイアモンド
		// const double nt = sqrt(1 + 4.3356 * lambda2/(lambda2 - pow(0.1060, 2)) + 0.3306 * lambda2/(lambda2 - pow(0.1750, 2)));

		const double nnt = into ? nc / nt : nt / nc;
		const double ddn = Dot(ray.dir, orienting_normal);
		const double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);

		if (cos2t < 0.0) { // 全反射した
			return (direct_light + obj.color * (radiance(reflection_ray, depth+1, wavelength))) / russian_roulette_probability;
		}
		// 屈折していく方向
		Vec tdir = Normalize(ray.dir * nnt - normal * (into ? 1.0 : -1.0) * (ddn * nnt + sqrt(cos2t)));

		// SchlickによるFresnelの反射係数の近似
		const double a = nt - nc, b = nt + nc;
		const double R0 = (a * a) / (b * b);
		const double c = 1.0 - (into ? -ddn : Dot(tdir, normal));
		const double Re = R0 + (1.0 - R0) * pow(c, 5.0);
		const double Tr = 1.0 - Re; // 屈折光の運ぶ光の量
		const double probability = 0.25 + 0.5 * Re;

		Ray refraction_ray = Ray(hitpoint, tdir);

		intersect_scene(refraction_ray, &lt, &lid);
		Spectrum direct_light_refraction(0.0);
		if (lid == LightID) {
			direct_light_refraction = light_intensity(wavelength);
		}

		// 一定以上レイを追跡したら屈折と反射のどちらか一方を追跡する。（さもないと指数的にレイが増える）
		// ロシアンルーレットで決定する。
		if (depth > MaxDepth) {
			if (rand01() < probability) { // 反射
				return 
					obj.color * (direct_light + radiance(reflection_ray, depth+1, wavelength) * Re)
					/ probability
					/ russian_roulette_probability;
			} else { // 屈折
				return 
					obj.color * (direct_light_refraction + radiance(refraction_ray, depth+1, wavelength)) * Tr
					/ (1.0 - probability)
					/ russian_roulette_probability;
			}
		} else { // 屈折と反射の両方を追跡
			return 
				obj.color * ((direct_light + radiance(reflection_ray, depth+1, wavelength)) * Re
				+ (direct_light_refraction + radiance(refraction_ray, depth+1, wavelength)) * Tr) / russian_roulette_probability;
		}
	} break;
	}

	return Spectrum(0.0);
}

// *** .hdrフォーマットで出力するための関数 ***
struct HDRPixel {
	unsigned char r, g, b, e;
	HDRPixel(const unsigned char r_ = 0, const unsigned char g_ = 0, const unsigned char b_ = 0, const unsigned char e_ = 0) :
	r(r_), g(g_), b(b_), e(e_) {};
	unsigned char get(int idx) {
		switch (idx) {
		case 0: return r;
		case 1: return g;
		case 2: return b;
		case 3: return e;
		} return 0;
	}

};

// doubleのRGB要素を.hdrフォーマット用に変換
HDRPixel get_hdr_pixel(const Color &color) {
	double d = std::max(color.x, std::max(color.y, color.z));
	if (d <= 1e-32)
		return HDRPixel();
	int e;
	double m = frexp(d, &e); // d = m * 2^e
	d = m * 256.0 / d;
	return HDRPixel(color.x * d, color.y * d, color.z * d, e + 128);
}

// 書き出し用関数
void save_hdr_file(const std::string &filename, const Color* image, const int width, const int height) {
	FILE *fp = fopen(filename.c_str(), "wb");
	if (fp == NULL) {
		std::cerr << "Error: " << filename << std::endl;
		return;
	}
	// .hdrフォーマットに従ってデータを書きだす
	// ヘッダ
	unsigned char ret = 0x0a;
	fprintf(fp, "#?RADIANCE%c", (unsigned char)ret);
	fprintf(fp, "# Made with 100%% pure HDR Shop%c", ret);
	fprintf(fp, "FORMAT=32-bit_rle_rgbe%c", ret);
	fprintf(fp, "EXPOSURE=1.0000000000000%c%c", ret, ret);

	// 輝度値書き出し
	fprintf(fp, "-Y %d +X %d%c", height, width, ret);
	for (int i = height - 1; i >= 0; i --) {
		std::vector<HDRPixel> line;
		for (int j = 0; j < width; j ++) {
			HDRPixel p = get_hdr_pixel(image[j + i * width]);
			line.push_back(p);
		}
		fprintf(fp, "%c%c", 0x02, 0x02);
		fprintf(fp, "%c%c", (width >> 8) & 0xFF, width & 0xFF);
		for (int i = 0; i < 4; i ++) {
			for (int cursor = 0; cursor < width;) {
				const int cursor_move = std::min(127, width - cursor);
				fprintf(fp, "%c", cursor_move);
				for (int j = cursor;  j < cursor + cursor_move; j ++)
					fprintf(fp, "%c", line[j].get(i));
				cursor += cursor_move;
			}
		}
	}

	fclose(fp);
}

// http://cvrl.ioo.ucl.ac.uk/cmfs.htm CIE 1964 10-deg, XYZ CMFs より
// 360nm - 830nm (5nm幅) で、波長->XYZ の関数
static const double wavelength2xyz_table[] = {
	1.222E-07,1.3398E-08,5.35027E-07, 
	9.1927E-07,1.0065E-07,4.0283E-06, 
	5.9586E-06,6.511E-07,2.61437E-05, 
	0.000033266,0.000003625,0.00014622, 
	0.000159952,0.000017364,0.000704776, 
	0.00066244,0.00007156,0.0029278, 
	0.0023616,0.0002534,0.0104822, 
	0.0072423,0.0007685,0.032344, 
	0.0191097,0.0020044,0.0860109, 
	0.0434,0.004509,0.19712, 
	0.084736,0.008756,0.389366, 
	0.140638,0.014456,0.65676, 
	0.204492,0.021391,0.972542, 
	0.264737,0.029497,1.2825, 
	0.314679,0.038676,1.55348, 
	0.357719,0.049602,1.7985, 
	0.383734,0.062077,1.96728, 
	0.386726,0.074704,2.0273, 
	0.370702,0.089456,1.9948, 
	0.342957,0.106256,1.9007, 
	0.302273,0.128201,1.74537, 
	0.254085,0.152761,1.5549, 
	0.195618,0.18519,1.31756, 
	0.132349,0.21994,1.0302, 
	0.080507,0.253589,0.772125, 
	0.041072,0.297665,0.57006, 
	0.016172,0.339133,0.415254, 
	0.005132,0.395379,0.302356, 
	0.003816,0.460777,0.218502, 
	0.015444,0.53136,0.159249, 
	0.037465,0.606741,0.112044, 
	0.071358,0.68566,0.082248, 
	0.117749,0.761757,0.060709, 
	0.172953,0.82333,0.04305, 
	0.236491,0.875211,0.030451, 
	0.304213,0.92381,0.020584, 
	0.376772,0.961988,0.013676, 
	0.451584,0.9822,0.007918, 
	0.529826,0.991761,0.003988, 
	0.616053,0.99911,0.001091, 
	0.705224,0.99734,0, 
	0.793832,0.98238,0, 
	0.878655,0.955552,0, 
	0.951162,0.915175,0, 
	1.01416,0.868934,0, 
	1.0743,0.825623,0, 
	1.11852,0.777405,0, 
	1.1343,0.720353,0, 
	1.12399,0.658341,0, 
	1.0891,0.593878,0, 
	1.03048,0.527963,0, 
	0.95074,0.461834,0, 
	0.856297,0.398057,0, 
	0.75493,0.339554,0, 
	0.647467,0.283493,0, 
	0.53511,0.228254,0, 
	0.431567,0.179828,0, 
	0.34369,0.140211,0, 
	0.268329,0.107633,0, 
	0.2043,0.081187,0, 
	0.152568,0.060281,0, 
	0.11221,0.044096,0, 
	0.0812606,0.0318004,0, 
	0.05793,0.0226017,0, 
	0.0408508,0.0159051,0, 
	0.028623,0.0111303,0, 
	0.0199413,0.0077488,0, 
	0.013842,0.0053751,0, 
	0.00957688,0.00371774,0, 
	0.0066052,0.00256456,0, 
	0.00455263,0.00176847,0, 
	0.0031447,0.00122239,0, 
	0.00217496,0.00084619,0, 
	0.0015057,0.00058644,0, 
	0.00104476,0.00040741,0, 
	0.00072745,0.000284041,0, 
	0.000508258,0.00019873,0, 
	0.00035638,0.00013955,0, 
	0.000250969,0.000098428,0, 
	0.00017773,0.000069819,0, 
	0.00012639,0.000049737,0, 
	0.000090151,3.55405E-05,0, 
	6.45258E-05,0.000025486,0, 
	0.000046339,1.83384E-05,0, 
	3.34117E-05,0.000013249,0, 
	0.000024209,9.6196E-06,0, 
	1.76115E-05,7.0128E-06,0, 
	0.000012855,5.1298E-06,0, 
	9.41363E-06,3.76473E-06,0, 
	0.000006913,2.77081E-06,0, 
	5.09347E-06,2.04613E-06,0, 
	3.7671E-06,1.51677E-06,0, 
	2.79531E-06,1.12809E-06,0, 
	0.000002082,8.4216E-07,0, 
	1.55314E-06,6.297E-07,0, 
};

// XYZ値をRGB値に変換する
// CIE RGB
Color xyz2rgb(const Color &xyz) {
	return Color(Dot(Color( 2.3655, -0.8971, -0.4683), xyz),
				 Dot(Color(-0.5151,  1.4264,  0.0887), xyz),
				 Dot(Color( 0.0052, -0.0144,  1.0089), xyz));
}

const int wl_max = 95;

int main(int argc, char **argv) {
	int width = 640;
	int height = 480;
	int samples = 32;

	// カメラ位置
	Ray camera(Vec(50.0, 52.0, 295.6), Normalize(Vec(0.0, -0.042612, -1.0)));
	// シーン内でのスクリーンのx,y方向のベクトル
	Vec cx = Vec(width * 0.5135 / height);
	Vec cy = Normalize(Cross(cx, camera.dir)) * 0.5135;
	Color *image = new Color[width * height];


	// 下で各波長についてパストレーシングするわけだが、波長の選択を完全にランダムにするのは効率が悪い。
	// 波長ごとに最終的な画像への寄与度は異なるからである。そこで、各波長の寄与度に基づいた重点サンプリングをする。
	// 具体的には各波長について、XYZ応答のうちY、すなわち輝度分に比例する確率密度関数を作り、それに基づいてサンプリングするようにする。
	// pdfはその密度関数でcdfはそこからサンプリングするための累積分布関数
	double cdf[wl_max];
	double pdf[wl_max];
	double luminance_table[wl_max];
	double total = 0.0;
	
	for (int i = 0; i < wl_max; i ++) {
		luminance_table[i] = wavelength2xyz_table[i * 3 + 1];
	}
		
	for (int i = 0; i < wl_max; i ++) {
		total += luminance_table[i];
		cdf[i] = total;
		pdf[i] = luminance_table[i];
	}
	for (int i = 0; i < wl_max; i ++) {
		cdf[i] /= total;
		pdf[i] /= total;
	}

// #pragma omp parallel for schedule(dynamic, 1) num_threads(3)
	for (int y = 0; y < height; y ++) {
		std::cerr << "Rendering (" << samples * 4 << " spp) " << (100.0 * y / (height - 1)) << "%" << std::endl;
		srand(y * y * y);
		for (int x = 0; x < width; x ++) {
			int image_index = y * width + x;	
			image[image_index] = Color();

			// 2x2のサブピクセルサンプリング
			for (int sy = 0; sy < 2; sy ++) {
				for (int sx = 0; sx < 2; sx ++) {
					Color accumulated_radiance = Color();
					// 一つのサブピクセルあたりsamples回サンプリングする
					for (int s = 0; s < samples; s ++) {
						// テントフィルターによってサンプリング
						// ピクセル範囲で一様にサンプリングするのではなく、ピクセル中央付近にサンプルがたくさん集まるように偏りを生じさせる
						const double r1 = 2.0 * rand01(), dx = r1 < 1.0 ? sqrt(r1) - 1.0 : 1.0 - sqrt(2.0 - r1);
						const double r2 = 2.0 * rand01(), dy = r2 < 1.0 ? sqrt(r2) - 1.0 : 1.0 - sqrt(2.0 - r2);
						Vec dir = cx * (((sx + 0.5 + dx) / 2.0 + x) / width - 0.5) +
								  cy * (((sy + 0.5 + dy) / 2.0 + y) / height- 0.5) + camera.dir;
						 
						// 波長を重点サンプリング
						double* p = std::lower_bound(cdf, cdf + wl_max, rand01());
						int wavelength_index = p - cdf;
						if (wavelength_index >= wl_max) wavelength_index = wl_max - 1;
						const double div = pdf[wavelength_index];
						const double wavelength = (wavelength_index * 5 + 360) / 1000.0;
						
						Spectrum value = radiance(Ray(camera.org + dir * 130.0, Normalize(dir)), 0, wavelength) / samples / div;
					
						accumulated_radiance = accumulated_radiance + xyz2rgb(
							Color(wavelength2xyz_table[wavelength_index * 3], 
							      wavelength2xyz_table[wavelength_index * 3 + 1],
							      wavelength2xyz_table[wavelength_index * 3 + 2]) * value);
					
					
					}
					image[image_index] = image[image_index] + accumulated_radiance;
				}
			}
		}
	}
	
	// .hdrフォーマットで出力
	save_hdr_file(std::string("image.hdr"), image, width, height);
}
