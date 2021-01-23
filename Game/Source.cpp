#define OLC_PGE_APPLICATION
#define _USE_MATH_DEFINES
#include "olcPixelGameEngine.h"


template<typename T>
struct v3d_generic {
	T x, y, z;

	v3d_generic() : x(0), y(0), z(0) {}
	v3d_generic(T _x, T _y, T _z) : x(_x), y(_y), z(_z) {}
	v3d_generic(const v3d_generic& v) : x(v.x), y(v.y), z(v.z) {}
	v3d_generic& operator=(const v3d_generic& v) = default;
	T mag() const { return T(std::sqrt(x * x + y * y + z * z)); }
	T mag2() const { return x * x + y * y + z * z; }
	v3d_generic  norm() const { T r = 1 / mag(); return v3d_generic(x * r, y * r, z * r); }
	//v2d_generic  perp() const { return v2d_generic(-y, x); }
	v3d_generic  floor() const { return v3d_generic(std::floor(x), std::floor(y), std::floor(z)); }
	v3d_generic  ceil() const { return v3d_generic(std::ceil(x), std::ceil(y), std::ceil(z)); }
	v3d_generic  max(const v3d_generic& v) const { return v3d_generic(std::max(x, v.x), std::max(y, v.y), std::max(z, v.z)); }
	v3d_generic  min(const v3d_generic& v) const { return v3d_generic(std::min(x, v.x), std::min(y, v.y), std::min(z, v.z)); }
	//T dot(const v2d_generic& rhs) const { return this->x * rhs.x + this->y * rhs.y; }
	//T cross(const v2d_generic& rhs) const { return this->x * rhs.y - this->y * rhs.x; }
	v3d_generic  operator +  (const v3d_generic& rhs) const { return v3d_generic(this->x + rhs.x, this->y + rhs.y, this->z + rhs.z); }
	v3d_generic  operator -  (const v3d_generic& rhs) const { return v3d_generic(this->x - rhs.x, this->y - rhs.y, this->z - rhs.z); }
	v3d_generic  operator *  (const T& rhs)           const { return v3d_generic(this->x * rhs, this->y * rhs, this->z * rhs); }
	v3d_generic  operator *  (const v3d_generic& rhs) const { return v3d_generic(this->x * rhs.x, this->y * rhs.y, this->z * rhs.z); }
	v3d_generic  operator /  (const T& rhs)           const { return v3d_generic(this->x / rhs, this->y / rhs, this->z / z); }
	v3d_generic  operator /  (const v3d_generic& rhs) const { return v3d_generic(this->x / rhs.x, this->y / rhs.y, this->z / rhs.z); }
	v3d_generic& operator += (const v3d_generic& rhs) { this->x += rhs.x; this->y += rhs.y; this->z += rhs.z; return *this; }
	v3d_generic& operator -= (const v3d_generic& rhs) { this->x -= rhs.x; this->y -= rhs.y; this->z -= rhs.z; return *this; }
	v3d_generic& operator *= (const T& rhs) { this->x *= rhs; this->y *= rhs; this->z *= rhs; return *this; }
	v3d_generic& operator /= (const T& rhs) { this->x /= rhs; this->y /= rhs; this->z /= rhs; return *this; }
	v3d_generic& operator *= (const v3d_generic& rhs) { this->x *= rhs.x; this->y *= rhs.y; this->z *= rhs.z; return *this; }
	v3d_generic& operator /= (const v3d_generic& rhs) { this->x /= rhs.x; this->y /= rhs.y; this->z /= rhs.z; return *this; }
	v3d_generic  operator +  () const { return { +x, +y, +z }; }
	v3d_generic  operator -  () const { return { -x, -y, -z }; }
	bool operator == (const v3d_generic& rhs) const { return (this->x == rhs.x && this->y == rhs.y && this->z == rhs.z); }
	bool operator != (const v3d_generic& rhs) const { return (this->x != rhs.x || this->y != rhs.y || this->z != rhs.z); }
	const std::string str() const { return std::string("(") + std::to_string(this->x) + "," + std::to_string(this->y) + "," + std::to_string(this->z) + ")"; }
	friend std::ostream& operator << (std::ostream& os, const v3d_generic& rhs) { os << rhs.str(); return os; }
	operator v3d_generic<int32_t>() const { return { static_cast<int32_t>(this->x), static_cast<int32_t>(this->y), static_cast<int32_t>(this->z) }; }
	operator v3d_generic<float>() const { return { static_cast<float>(this->x), static_cast<float>(this->y), static_cast<float>(this->z) }; }
	operator v3d_generic<double>() const { return { static_cast<double>(this->x), static_cast<double>(this->y), static_cast<double>(this->z) }; }
};

typedef v3d_generic<int32_t> vi3d;
typedef v3d_generic<uint32_t> vu3d;
typedef v3d_generic<float> vf3d;
typedef v3d_generic<double> vd3d;


const uint32_t emptyVoxel = 0;


struct svo_base {

};


struct svo_object {
	uint8_t size;	// ?
	olc::Pixel color = olc::Pixel(emptyVoxel);
	std::vector<svo_object> childs;
	svo_object* parent;
	uint8_t used_childs_mask;

	svo_object* create_childs() {
		childs.resize(8);
		childs[0].parent = this;
		childs[1].parent = this;
		childs[2].parent = this;
		childs[3].parent = this;
		childs[4].parent = this;
		childs[5].parent = this;
		childs[6].parent = this;
		childs[7].parent = this;
	}
};


class game_object {
	vu3d size;
	vd3d pos;
};


template<typename T>
struct octree_node {
	olc::Pixel color = olc::Pixel(emptyVoxel);
	std::vector<T> items;
	std::vector<octree_node> childs;
	octree_node* parent = nullptr;
	uint32_t depth = 1;
	uint8_t used_childs_mask;

	bool create_childs() {
		if (depth == 1)
			return false;

		uint32_t child_depth = depth - 1;
		childs.clear();
		childs.resize(8);
		childs[0].parent = this;
		childs[0].depth = child_depth;
		childs[1].parent = this;
		childs[1].depth = child_depth;
		childs[2].parent = this;
		childs[2].depth = child_depth;
		childs[3].parent = this;
		childs[3].depth = child_depth;
		childs[4].parent = this;
		childs[4].depth = child_depth;
		childs[5].parent = this;
		childs[5].depth = child_depth;
		childs[6].parent = this;
		childs[6].depth = child_depth;
		childs[7].parent = this;
		childs[7].depth = child_depth;

		return true;
	}
	//void create_depth(uint32_t depth) {
	//	if (depth == 0)
	//		return;
	//	depth -= 1;
	//	create_childs();
	//	childs[0].create_depth(depth);
	//	childs[1].create_depth(depth);
	//	childs[2].create_depth(depth);
	//	childs[3].create_depth(depth);
	//	childs[4].create_depth(depth);
	//	childs[5].create_depth(depth);
	//	childs[6].create_depth(depth);
	//	childs[7].create_depth(depth);
	//}
};


class Space3D {
	octree_node<game_object> space;
	vu3d size;
	vi3d pos;

public:
	Space3D(vi3d&& coords, vu3d&& sizes) {
		size = sizes;
		pos = coords;

		uint32_t maxSide = std::max(size.x, std::max(size.y, size.z));
		uint32_t depth = 1;
		uint32_t boxSide = 1;
		while (true) {
			if (boxSide >= maxSide)
				break;
			depth += 1;
			boxSide <<= 1;
		}
		space.depth = depth;
	}
};


class camera {
	double x, y, z;
	double hAngle, vAngle;
	const double move_speed = 0.1;
	const double turn_speed = 0.001;
	const double fov = M_PI / 2;	// 90 degree

public:
	camera(double _x, double _y, double _z)
		: x(_x), y(_y), z(_z), hAngle(0.0), vAngle(0.0) {}

public:
	void turnRight() {
		hAngle += turn_speed;
		if (hAngle > M_PI) hAngle -= 2.0 * M_PI;
	}
	void turnLeft() {
		hAngle -= turn_speed;
		if (hAngle < -M_PI) hAngle += 2.0 * M_PI;
	}
	void turnUp() {
		if (vAngle <= M_PI / 2.0)
			vAngle += turn_speed;
	}
	void turnDown() {
		if (vAngle >= -M_PI / 2.0)
			vAngle -= turn_speed;
	}

	void moveW() {
		x += move_speed * std::cos(hAngle) * std::cos(vAngle);
		y += move_speed                    * std::sin(vAngle);
		z += move_speed * std::sin(hAngle) * std::cos(vAngle);
	}
	void moveA() {
		double _hAngle = hAngle + M_PI / 2.0;
		x += move_speed * std::cos(_hAngle) * std::cos(vAngle);
		y += move_speed * std::sin(vAngle);
		z += move_speed * std::sin(_hAngle) * std::cos(vAngle);
	}
	void moveS() {
		x -= move_speed * std::cos(hAngle) * std::cos(vAngle);
		y -= move_speed                    * std::sin(vAngle);
		z -= move_speed * std::sin(hAngle) * std::cos(vAngle);
	}
	void moveD() {
		double _hAngle = hAngle - M_PI / 2.0;
		x += move_speed * std::cos(_hAngle) * std::cos(vAngle);
		y += move_speed * std::sin(vAngle);
		z += move_speed * std::sin(_hAngle) * std::cos(vAngle);
	}

	double X() const {
		return x;
	}
	double Y() const {
		return y;
	}
	double Z() const {
		return z;
	}

	vd3d norm_direction() const {
		return {
			std::cos(hAngle) * std::cos(vAngle),
			std::sin(vAngle),
			std::sin(hAngle) * std::cos(vAngle),
		};
	}
};


class SVOGame : public olc::PixelGameEngine {
public:
	SVOGame() {
		sAppName = "SVO";
	}


protected:
	//game_object test_object;
	camera player_view = camera(-10.0, 0.0, 0.0);
	Space3D sapce = { { 0, 0, 0 }, { 16, 16, 16 } };


public:
	bool OnUserCreate() override {
		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override {
		// Movement
		if (GetKey(olc::W).bHeld)
			player_view.moveW();
		if (GetKey(olc::A).bHeld)
			player_view.moveA();
		if (GetKey(olc::S).bHeld)
			player_view.moveS();
		if (GetKey(olc::D).bHeld)
			player_view.moveD();
		if (GetKey(olc::Q).bHeld)
			player_view.turnLeft();
		if (GetKey(olc::E).bHeld)
			player_view.turnRight();


		// Render
		Clear(olc::BLACK);

		/*
		for pixel:
			
		*/

		std::stringstream coords;
		coords << "X: " << player_view.X() << '\n';
		coords << "Y: " << player_view.Y() << '\n';
		coords << "Z: " << player_view.Z() << '\n';
		coords << "dir: " << player_view.norm_direction().str() << '\n';
		DrawString(0, 0, coords.str());


		return true;
	}


protected:
};


int main(int argc, char* argv[]) {
	SVOGame game;
	if (game.Construct(480, 360, 2, 2))
		game.Start();
	return 0;
}