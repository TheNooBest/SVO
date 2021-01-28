#define OLC_PGE_APPLICATION
#define _USE_MATH_DEFINES
#include "olcPixelGameEngine.h"
#include <stack>


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


struct svo_object {
	olc::Pixel color = olc::Pixel(emptyVoxel);
	std::vector<svo_object> childs;
	svo_object* parent = nullptr;
	uint32_t depth = 1;
	uint8_t used_childs_mask = 0x00;

	void set_used_child(svo_object* child) {
		if (child == &childs[0]) {
			used_childs_mask |= 0x01;
			return;
		}
		if (child == &childs[1]) {
			used_childs_mask |= 0x02;
			return;
		}
		if (child == &childs[2]) {
			used_childs_mask |= 0x04;
			return;
		}
		if (child == &childs[3]) {
			used_childs_mask |= 0x08;
			return;
		}
		if (child == &childs[4]) {
			used_childs_mask |= 0x10;
			return;
		}
		if (child == &childs[5]) {
			used_childs_mask |= 0x20;
			return;
		}
		if (child == &childs[6]) {
			used_childs_mask |= 0x40;
			return;
		}
		if (child == &childs[7]) {
			used_childs_mask |= 0x80;
			return;
		}
	}
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

		if (parent != nullptr)
			parent->set_used_child(this);

		return true;
	}
};


class GameObject {
	vu3d size;
	vd3d pos;
};


template<typename T>
struct octree_node {
	std::vector<T> items;
	svo_object svo;
};


class Space3D {
	octree_node<GameObject> space;
	vu3d size;
	vi3d pos;
	vi3d oppositePos;

public:
	Space3D() {}
	Space3D(vi3d&& coords, vu3d&& sizes) {
		size = sizes;
		pos = coords;
		oppositePos = pos + size;

		uint32_t maxSide = std::max(size.x, std::max(size.y, size.z));
		uint32_t depth = 1;
		uint32_t boxSide = 1;
		while (boxSide < maxSide) {
			depth += 1;
			boxSide <<= 1;
		}
		space.svo.depth = depth;
	}

	void InitFrom(svo_object& svo) {
		// todo: sizes
		space.svo = svo;
	}

	vu3d Size() const { return size; }
	vi3d Pos() const { return pos; }
	vi3d OPos() const { return oppositePos; }

	uint32_t Depth() const { return space.svo.depth; }
};


class Camera {

public:
	vd3d pos;
	double hAngle, vAngle;
	const double move_speed = 0.1;
	const double turn_speed = 0.001;
	const double fovH = M_PI / 2;	// 90 degree
	const double fovV = M_PI / 3;	// 60 degree

public:
	Camera()
		: pos(0.0, 0.0, 0.0), hAngle(0.0), vAngle(0.0) {}
	Camera(double _x, double _y, double _z)
		: pos(_x, _y, _z), hAngle(0.0), vAngle(0.0) {}

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
		pos.x += move_speed * std::cos(hAngle) * std::cos(vAngle);
		pos.y += move_speed                    * std::sin(vAngle);
		pos.z += move_speed * std::sin(hAngle) * std::cos(vAngle);
	}
	void moveA() {
		double _hAngle = hAngle + M_PI / 2.0;
		pos.x += move_speed * std::cos(_hAngle) * std::cos(vAngle);
		pos.y += move_speed * std::sin(vAngle);
		pos.z += move_speed * std::sin(_hAngle) * std::cos(vAngle);
	}
	void moveS() {
		pos.x -= move_speed * std::cos(hAngle) * std::cos(vAngle);
		pos.y -= move_speed                    * std::sin(vAngle);
		pos.z -= move_speed * std::sin(hAngle) * std::cos(vAngle);
	}
	void moveD() {
		double _hAngle = hAngle - M_PI / 2.0;
		pos.x += move_speed * std::cos(_hAngle) * std::cos(vAngle);
		pos.y += move_speed * std::sin(vAngle);
		pos.z += move_speed * std::sin(_hAngle) * std::cos(vAngle);
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
	Camera player_view = Camera(-10.0, 0.0, 0.0);
	Space3D space;
	std::vector<double> vAngles;
	std::vector<olc::vd2d> xzComponents;
	std::vector<octree_node<GameObject>*> octree_node_stack;
	uint32_t st_index;


protected:
	bool OnUserCreate() override {
		svo_object svo;
		svo.depth = 2;
		svo.create_childs();
		svo.childs[0].color = olc::RED;
		svo.childs[2].color = olc::RED;
		space.InitFrom(svo);

		vAngles.resize(ScreenHeight());
		xzComponents.resize(ScreenWidth());

		octree_node_stack.resize(space.Depth() * (size_t)4);

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

		double leftAngle = player_view.hAngle - player_view.fovH / 2.0;
		double rightAngle = player_view.hAngle + player_view.fovH / 2.0;
		double bottomAngle = player_view.vAngle - player_view.fovV / 2.0;
		double upperAngle = player_view.vAngle + player_view.fovV / 2.0;

		double vDelta = (upperAngle - bottomAngle) / ScreenHeight();
		double hDelta = (rightAngle - leftAngle) / ScreenWidth();
		for (int i = 0; i < ScreenHeight(); i++) {
			vAngles[i] = bottomAngle + i * vDelta;
		}
		for (int i = 0; i < ScreenWidth(); i++) {
			double hAngle = leftAngle + i * hDelta;
			xzComponents[i] = { cos(hAngle), sin(hAngle) };
		}

		for (double vAngle : vAngles) {
			double rayY = sin(vAngle);
			double xzLength = cos(vAngle);
			for (olc::vd2d& xzDir : xzComponents) {
				vd3d rayDir = { xzDir.x * xzLength, rayY, xzDir.y * xzLength };
				//                                     ^^^^^^^ <- zComponent
				vd3d raySource = player_view.pos + rayDir;

				calcRay(raySource, rayDir);
			}
		}

		std::stringstream coords;
		coords << "pos: " << player_view.pos.str() << '\n';
		coords << "dir: " << player_view.norm_direction().str() << '\n';
		DrawString(0, 0, coords.str());


		return true;
	}


protected:
	inline void calcRay(vd3d& raySource, vd3d& rayDir) {
		vd3d invRay = vd3d{ 1, 1, 1 } / rayDir;

		double t1 = (space.Pos().x - raySource.x) * invRay.x;
		double t2 = (space.OPos().x - raySource.x) * invRay.x;
		double tmin = std::min(t1, t2);
		double tmax = std::max(t1, t2);

		t1 = (space.Pos().y - raySource.y) * invRay.y;
		t2 = (space.OPos().y - raySource.y) * invRay.y;
		tmin = std::max(tmin, std::min(t1, t2));
		tmax = std::min(tmax, std::max(t1, t2));

		t1 = (space.Pos().z - raySource.z) * invRay.z;
		t2 = (space.OPos().z - raySource.z) * invRay.z;
		tmin = std::max(tmin, std::min(t1, t2));
		tmax = std::min(tmax, std::max(t1, t2));

		if (tmin > tmax)
			return;

		// 1. find collisions with global octree
		// 2.1. if empty -> next ray
		// 2.2. put in stack (first node in path of ray on top of stack, last on bottom)
		// 3. check top node
		// 3.1. if empty -> next node
		// 3.1.1. if stack is empty -> next ray
		// 3.2. if collide -> go deeper or render (may depends on distance)
		// ?. rays to light sources
		// ?.?. reflections...
		// 4. next ray
	}


protected:
};


int main(int argc, char* argv[]) {
	SVOGame game;
	if (game.Construct(480, 360, 2, 2))
		game.Start();
	return 0;
}