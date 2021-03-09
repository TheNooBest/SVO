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
	T dot(const v3d_generic& rhs) const { return this->x * rhs.x + this->y * rhs.y + this->z * rhs.z; }
	v3d_generic  cross(const v3d_generic& rhs) const { return { this->y * rhs.z - this->z * rhs.y, this->z * rhs.x - this->x * rhs.z, this->x * rhs.y - this->y * rhs.z }; }
	v3d_generic  operator +  (const T& rhs)			  const { return v3d_generic(this->x + rhs, this->y + rhs, this->z + rhs); }
	v3d_generic  operator +  (const v3d_generic& rhs) const { return v3d_generic(this->x + rhs.x, this->y + rhs.y, this->z + rhs.z); }
	v3d_generic  operator -  (const T& rhs)			  const { return v3d_generic(this->x - rhs, this->y - rhs, this->z - rhs); }
	v3d_generic  operator -  (const v3d_generic& rhs) const { return v3d_generic(this->x - rhs.x, this->y - rhs.y, this->z - rhs.z); }
	v3d_generic  operator *  (const T& rhs)           const { return v3d_generic(this->x * rhs, this->y * rhs, this->z * rhs); }
	v3d_generic  operator *  (const v3d_generic& rhs) const { return v3d_generic(this->x * rhs.x, this->y * rhs.y, this->z * rhs.z); }
	v3d_generic  operator /  (const T& rhs)           const { return v3d_generic(this->x / rhs, this->y / rhs, this->z / rhs); }
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

	//custom
	v3d_generic sum() const { return this->x + this->y + this->z; }
	v3d_generic rotate(v3d_generic& a, double angle) {
		double c = std::cos(angle), s = std::sin(angle), cc = 1 - c;
		double xy = a.x * a.y, xz = a.x * a.z, yz = a.y * a.z;
		return {
			(c + cc * a.x * a.x) * x + (cc * xy - s * a.z)  * y + (cc * xz + s * a.y)  * z,
			(cc * xy + s * a.z)  * x + (c + cc * a.y * a.y) * y + (cc * yz - s * a.x)  * z,
			(cc * xz - s * a.y)  * x + (cc * yz + s * a.x)  * y + (c + cc * a.z * a.z) * z
		};
		//v3d_generic<T> k = a.norm();
		//double c = std::cos(angle), s = std::sin(angle);
		//return (*this * c) + k.cross(*this) * s + k * k.dot(*this) * (1 - c);
	}
};

typedef v3d_generic<int32_t> vi3d;
typedef v3d_generic<uint32_t> vu3d;
typedef v3d_generic<float> vf3d;
typedef v3d_generic<double> vd3d;


class GameObject;


struct svo_space {
	std::vector<GameObject> items;
	std::vector<svo_space> childs;
	svo_space* parent = nullptr;
	uint32_t depth = 1;
	uint8_t used_childs_mask = 0x00;

	void set_used_child(svo_space* child) {
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


struct svo_model {
	olc::Pixel color = olc::BLANK;
	std::vector<svo_model> childs;
	svo_model* parent = nullptr;
	uint32_t depth = 1;
	uint8_t used_childs_mask = 0x00;

	void set_used_child(svo_model* child) {
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

	bool is_terminal() {
		return childs.empty();
	}
};


class GameObject {
public:
	vu3d size;
	vd3d pos;
	svo_model model;
};


class Space3D {
public:
	svo_space svo;
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
		svo.depth = depth;
	}


	vu3d Size() const { return size; }
	vi3d Pos() const { return pos; }
	vi3d OPos() const { return oppositePos; }

	uint32_t Depth() const { return svo.depth; }
};


class Camera {

public:
	vd3d pos;
	double hAngle, vAngle;
	const double move_speed = 0.1;
	const double turn_speed = 0.01;
	const double fovH = M_PI / 2;	// 90 degree
	const double fovV = M_PI / 3;	// 60 degree

public:
	Camera()
		: pos(0.0, 0.0, 0.0), hAngle(0.0), vAngle(0.0) {}
	Camera(double _x, double _y, double _z)
		: pos(_x, _y, _z), hAngle(0.0), vAngle(0.0) {}
	Camera(double _x, double _y, double _z, double _hAngle, double _vAngle)
		: pos(_x, _y, _z), hAngle(_hAngle), vAngle(_vAngle) {}

public:
	void turnRight() {
		hAngle -= turn_speed;
		if (hAngle < -M_PI) hAngle += M_PI + M_PI;
	}
	void turnLeft() {
		hAngle += turn_speed;
		if (hAngle > M_PI) hAngle -= M_PI + M_PI;
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
		pos.z += move_speed * -std::sin(hAngle) * std::cos(vAngle);
	}
	void moveA() {
		double _hAngle = hAngle + M_PI / 2.0;
		pos.x += move_speed * std::cos(_hAngle) * std::cos(vAngle);
		pos.y += move_speed * std::sin(vAngle);
		pos.z += move_speed * -std::sin(_hAngle) * std::cos(vAngle);
	}
	void moveS() {
		pos.x -= move_speed * std::cos(hAngle) * std::cos(vAngle);
		pos.y -= move_speed                    * std::sin(vAngle);
		pos.z -= move_speed * -std::sin(hAngle) * std::cos(vAngle);
	}
	void moveD() {
		double _hAngle = hAngle - M_PI / 2.0;
		pos.x += move_speed * std::cos(_hAngle) * std::cos(vAngle);
		pos.y += move_speed * std::sin(vAngle);
		pos.z += move_speed * -std::sin(_hAngle) * std::cos(vAngle);
	}
	void moveUp() {
		pos.x += move_speed * -std::sin(vAngle) * -std::cos(hAngle);
		pos.y += move_speed * std::cos(vAngle);
		pos.z += move_speed * -std::sin(vAngle) * -std::sin(hAngle);
	}
	void moveDown() {
		pos.x -= move_speed * -std::sin(vAngle) * -std::cos(hAngle);
		pos.y -= move_speed * std::cos(vAngle);
		pos.z -= move_speed * -std::sin(vAngle) * -std::sin(hAngle);
	}

	vd3d norm_direction() const {
		return {
			std::cos(hAngle) * std::cos(vAngle),
			std::sin(vAngle),
			-std::sin(hAngle) * std::cos(vAngle),
		};
	}
	vd3d perp_vert() const {
		double _vAngle = vAngle + M_PI / 2.0;
		return {
			std::cos(hAngle) * std::cos(_vAngle),
			std::sin(_vAngle),
			-std::sin(hAngle) * std::cos(_vAngle),
		};
	}
	vd3d perp_hor() const {
		double _hAngle = hAngle + M_PI / 2.0;
		return {
			-std::cos(_hAngle) * std::cos(vAngle),
			std::sin(vAngle),
			std::sin(_hAngle) * std::cos(vAngle),
		};
	}
};


class SVOGame : public olc::PixelGameEngine {
public:
	SVOGame() {
		sAppName = "SVO";
	}


protected:
	Camera player_view = Camera(-1.0, 0.0, 0.0);
	Space3D space;
	GameObject testObject;
	uint32_t st_index;
	std::vector<double> vAngles;
	std::vector<double> hAngles;
	std::vector<vd3d>   horAxis;


protected:
	bool OnUserCreate() override {
		svo_model model;
		model.depth = 3;
		model.create_childs();
		model.childs[0].create_childs();
		model.childs[1].color = olc::GREEN;
		model.childs[2].color = olc::YELLOW;
		model.childs[3].color = olc::BLUE;
		model.childs[4].color = olc::DARK_RED;
		model.childs[5].color = olc::DARK_GREEN;
		model.childs[6].color = olc::DARK_YELLOW;
		model.childs[7].color = olc::DARK_BLUE;

		model.childs[0].childs[0].color = olc::RED;
		model.childs[0].childs[1].color = olc::GREEN;
		model.childs[0].childs[2].color = olc::YELLOW;
		model.childs[0].childs[3].color = olc::BLUE;
		model.childs[0].childs[4].color = olc::DARK_RED;
		model.childs[0].childs[5].color = olc::DARK_GREEN;
		model.childs[0].childs[6].color = olc::DARK_YELLOW;
		model.childs[0].childs[7].color = olc::DARK_BLUE;

		testObject.model = model;
		testObject.pos = { 0, 0, 0 };
		uint32_t size = 1 << (model.depth - 1);
		testObject.size = { size, size, size };

		vAngles.resize(ScreenHeight());
		hAngles.resize(ScreenWidth());
		horAxis.resize(ScreenWidth());

		double rightAngle = -player_view.fovH / 2.0;
		double leftAngle = player_view.fovH / 2.0;
		double bottomAngle = -player_view.fovV / 2.0;
		double upperAngle = player_view.fovV / 2.0;

		double vDelta = (upperAngle - bottomAngle) / ScreenHeight();
		double hDelta = (leftAngle - rightAngle) / ScreenWidth();
		for (int i = 0; i < vAngles.size(); i++)
			vAngles[i] = upperAngle - i * vDelta;
		for (int i = 0; i < hAngles.size(); i++)
			hAngles[i] = leftAngle - i * hDelta;

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
		if (GetKey(olc::SPACE).bHeld)
			player_view.moveUp();
		if (GetKey(olc::SHIFT).bHeld)
			player_view.moveDown();
		if (GetKey(olc::Q).bHeld)
			player_view.turnLeft();
		if (GetKey(olc::E).bHeld)
			player_view.turnRight();


		// Render
		Clear(olc::BLACK);

		vd3d vert = player_view.perp_vert() * std::tan(player_view.fovV / 2.0);
		vd3d hor  = player_view.perp_hor()  * std::tan(player_view.fovH / 2.0);

		vd3d luCorner = player_view.pos + player_view.norm_direction() + vert - hor;
		vd3d ruCorner = player_view.pos + player_view.norm_direction() + vert + hor;
		vd3d lbCorner = player_view.pos + player_view.norm_direction() - vert - hor;
		vd3d rbCorner = player_view.pos + player_view.norm_direction() - vert + hor;

		vd3d deltaY = (lbCorner - luCorner) / ScreenHeight();
		vd3d deltaX = (ruCorner - luCorner) / ScreenWidth();

		for (int y = 0; y < ScreenHeight(); y++) {
			vd3d left = luCorner + deltaY * y;
			for (int x = 0; x < ScreenWidth(); x++) {
				vd3d raySource = left + deltaX * x;
				vd3d rayDir = raySource - player_view.pos;

				olc::Pixel pix = rayParameter(&testObject.model, testObject.pos, testObject.pos + testObject.size, raySource, rayDir);
				Draw(x, y, pix);
			}
		}

		olc::vi2d mousePos = GetMousePos();

		std::stringstream coords;
		coords << "pos: " << player_view.pos.str() << '\n';
		coords << "dir: " << player_view.norm_direction().str() << '\n';
		coords << "mouse dir: " << (luCorner + deltaY * mousePos.y + deltaX * mousePos.x - player_view.pos).str() << '\n';
		DrawString(0, 0, coords.str());


		return true;
	}


protected:
	olc::Pixel rayParameter(svo_model* node, vd3d node_0, vd3d node_1, vd3d ray_source, vd3d ray_dir) {
		unsigned char a = 0;

		ray_dir += node_0;
		node_1 -= node_0;

		uint32_t node_size = node_1.x;

		if (ray_dir.x < 0.0) {
			ray_source.x = node_size - ray_source.x;
			ray_dir.x = -ray_dir.x;
			a |= 1;
		}
		if (ray_dir.y < 0.0) {
			ray_source.y = node_size - ray_source.y;
			ray_dir.y = -ray_dir.y;
			a |= 2;
		}
		if (ray_dir.z < 0.0) {
			ray_source.z = node_size - ray_source.z;
			ray_dir.z = -ray_dir.z;
			a |= 4;
		}

		double divX = 1.0 / ray_dir.x;
		double divY = 1.0 / ray_dir.y;
		double divZ = 1.0 / ray_dir.z;

		double tx0 = -ray_source.x * divX;
		double tx1 = (node_1.x - ray_source.x) * divX;
		double ty0 = -ray_source.y * divY;
		double ty1 = (node_1.y - ray_source.y) * divY;
		double tz0 = -ray_source.z * divZ;
		double tz1 = (node_1.z - ray_source.z) * divZ;

		if (std::isnan(tx0))
			tx0 = -INFINITY;
		if (std::isnan(tx1))
			return olc::BLANK;
		if (std::isnan(ty0))
			ty0 = -INFINITY;
		if (std::isnan(ty1))
			return olc::BLANK;
		if (std::isnan(tz0))
			tz0 = -INFINITY;
		if (std::isnan(tz1))
			return olc::BLANK;

		if (std::max(tx0, std::max(ty0, tz0)) < std::min(tx1, std::min(ty1, tz1)))
			return procSubtree(tx0, ty0, tz0, tx1, ty1, tz1, node, a, node_1, ray_source);

		return olc::BLANK;
	}


	olc::Pixel procSubtree(double tx0, double ty0, double tz0, double tx1, double ty1, double tz1, svo_model* node, uint8_t a, vd3d node_1, vd3d ray_source) {
		double txm, tym, tzm;

		if (tx1 < 0.0 || ty1 < 0.0 || tz1 < 0.0)
			return olc::BLANK;
		if (node->is_terminal())
			return node->color;

		if (std::isinf(tx0))
			txm = ray_source.x < node_1.x * 0.5 ? INFINITY : -INFINITY;
		else
			txm = 0.5 * (tx0 + tx1);
		if (std::isinf(ty0))
			tym = ray_source.y < node_1.y * 0.5 ? INFINITY : -INFINITY;
		else
			tym = 0.5 * (ty0 + ty1);
		if (std::isinf(tz0))
			tzm = ray_source.z < node_1.z * 0.5 ? INFINITY : -INFINITY;
		else
			tzm = 0.5 * (tz0 + tz1);

		uint8_t currNode = firstNode(tx0, ty0, tz0, txm, tym, tzm);

		vd3d half_node = node_1 * 0.5;

		do {
			olc::Pixel pix;
			switch (currNode)
			{
			case 0:
				pix = procSubtree(tx0, ty0, tz0, txm, tym, tzm, &node->childs[a], a, half_node, ray_source);
				if (pix != olc::BLANK) return pix;
				currNode = nextNode(txm, 1, tym, 2, tzm, 4);
				break;
			case 1:
				pix = procSubtree(txm, ty0, tz0, tx1, tym, tzm, &node->childs[1 ^ a], a, half_node, ray_source + vd3d{ half_node.x, 0.0, 0.0 });
				if (pix != olc::BLANK) return pix;
				currNode = nextNode(tx1, 8, tym, 3, tzm, 5);
				break;
			case 2:
				pix = procSubtree(tx0, tym, tz0, txm, ty1, tzm, &node->childs[2 ^ a], a, half_node, ray_source + vd3d{ 0.0, half_node.y, 0.0 });
				if (pix != olc::BLANK) return pix;
				currNode = nextNode(txm, 3, ty1, 8, tzm, 6);
				break;
			case 3:
				pix = procSubtree(txm, tym, tz0, tx1, ty1, tzm, &node->childs[3 ^ a], a, half_node, ray_source + vd3d{ half_node.x, half_node.y, 0.0 });
				if (pix != olc::BLANK) return pix;
				currNode = nextNode(tx1, 8, ty1, 8, tzm, 7);
				break;
			case 4:
				pix = procSubtree(tx0, ty0, tzm, txm, tym, tz1, &node->childs[4 ^ a], a, half_node, ray_source + vd3d{ half_node.x, 0.0, half_node.z });
				if (pix != olc::BLANK) return pix;
				currNode = nextNode(txm, 5, tym, 6, tz1, 8);
				break;
			case 5:
				pix = procSubtree(txm, ty0, tzm, tx1, tym, tz1, &node->childs[5 ^ a], a, half_node, ray_source + vd3d{ half_node.x, 0.0, half_node.z });
				if (pix != olc::BLANK) return pix;
				currNode = nextNode(tx1, 8, tym, 7, tz1, 8);
				break;
			case 6:
				pix = procSubtree(tx0, tym, tzm, txm, ty1, tz1, &node->childs[6 ^ a], a, half_node, ray_source + vd3d{ half_node.x, half_node.y, half_node.z });
				if (pix != olc::BLANK) return pix;
				currNode = nextNode(txm, 7, ty1, 8, tz1, 8);
				break;
			case 7:
				pix = procSubtree(txm, tym, tzm, tx1, ty1, tz1, &node->childs[7 ^ a], a, half_node, ray_source + vd3d{ half_node.x, half_node.y, half_node.z });
				if (pix != olc::BLANK) return pix;
				currNode = 8;
				break;
			}
		} while (currNode < 8);

		return olc::BLANK;
	}


	uint8_t firstNode(double tx0, double ty0, double tz0, double txm, double tym, double tzm) {
		uint8_t answer = 0;
		if (tx0 > ty0) {
			if (tx0 > tz0) {
				if (tym < tx0) answer |= 2;
				if (tzm < tx0) answer |= 4;
				return answer;
			}
		}
		else {
			if (ty0 > tz0) {
				if (txm < ty0) answer |= 1;
				if (tzm < ty0) answer |= 4;
				return answer;
			}
		}
		if (txm < tz0) answer |= 1;
		if (tym < tz0) answer |= 2;
		return answer;
	}


	uint32_t nextNode(double txm, uint8_t x, double tym, uint8_t y, double tzm, uint8_t z) {
		if (txm < tym) {
			if (txm < tzm) return x;
		}
		else {
			if (tym < tzm) return y;
		}
		return z;
	}
};


int main(int argc, char* argv[]) {
	SVOGame game;
	if (game.Construct(600, 400, 2, 2))
		game.Start();
	return 0;
}