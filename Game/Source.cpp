#define OLC_PGE_APPLICATION
#define _USE_MATH_DEFINES
#define SVO_RENDER_V2
#include "olcPixelGameEngine.h"
#include <thread>
#include <condition_variable>


using namespace std::chrono_literals;


template <typename T>
struct v3d_generic
{
	T x, y, z;

	v3d_generic() : x(0), y(0), z(0) {}
	v3d_generic(T _x, T _y, T _z) : x(_x), y(_y), z(_z) {}
	v3d_generic(const v3d_generic &v) : x(v.x), y(v.y), z(v.z) {}
	v3d_generic &operator=(const v3d_generic &v) = default;
	T mag() const { return T(std::sqrt(x * x + y * y + z * z)); }
	T mag2() const { return x * x + y * y + z * z; }
	v3d_generic norm() const
	{
		T r = 1 / mag();
		return v3d_generic(x * r, y * r, z * r);
	}
	//v2d_generic  perp() const { return v2d_generic(-y, x); }
	v3d_generic floor() const { return v3d_generic(std::floor(x), std::floor(y), std::floor(z)); }
	v3d_generic ceil() const { return v3d_generic(std::ceil(x), std::ceil(y), std::ceil(z)); }
	v3d_generic max(const v3d_generic &v) const { return v3d_generic(std::max(x, v.x), std::max(y, v.y), std::max(z, v.z)); }
	v3d_generic min(const v3d_generic &v) const { return v3d_generic(std::min(x, v.x), std::min(y, v.y), std::min(z, v.z)); }
	T dot(const v3d_generic &rhs) const { return this->x * rhs.x + this->y * rhs.y + this->z * rhs.z; }
	v3d_generic cross(const v3d_generic &rhs) const { return {this->y * rhs.z - this->z * rhs.y, this->z * rhs.x - this->x * rhs.z, this->x * rhs.y - this->y * rhs.z}; }
	v3d_generic operator+(const T &rhs) const { return v3d_generic(this->x + rhs, this->y + rhs, this->z + rhs); }
	v3d_generic operator+(const v3d_generic &rhs) const { return v3d_generic(this->x + rhs.x, this->y + rhs.y, this->z + rhs.z); }
	v3d_generic operator-(const T &rhs) const { return v3d_generic(this->x - rhs, this->y - rhs, this->z - rhs); }
	v3d_generic operator-(const v3d_generic &rhs) const { return v3d_generic(this->x - rhs.x, this->y - rhs.y, this->z - rhs.z); }
	v3d_generic operator*(const T &rhs) const { return v3d_generic(this->x * rhs, this->y * rhs, this->z * rhs); }
	v3d_generic operator*(const v3d_generic &rhs) const { return v3d_generic(this->x * rhs.x, this->y * rhs.y, this->z * rhs.z); }
	v3d_generic operator/(const T &rhs) const { return v3d_generic(this->x / rhs, this->y / rhs, this->z / rhs); }
	v3d_generic operator/(const v3d_generic &rhs) const { return v3d_generic(this->x / rhs.x, this->y / rhs.y, this->z / rhs.z); }
	v3d_generic &operator+=(const v3d_generic &rhs)
	{
		this->x += rhs.x;
		this->y += rhs.y;
		this->z += rhs.z;
		return *this;
	}
	v3d_generic &operator-=(const v3d_generic &rhs)
	{
		this->x -= rhs.x;
		this->y -= rhs.y;
		this->z -= rhs.z;
		return *this;
	}
	v3d_generic &operator*=(const T &rhs)
	{
		this->x *= rhs;
		this->y *= rhs;
		this->z *= rhs;
		return *this;
	}
	v3d_generic &operator/=(const T &rhs)
	{
		this->x /= rhs;
		this->y /= rhs;
		this->z /= rhs;
		return *this;
	}
	v3d_generic &operator*=(const v3d_generic &rhs)
	{
		this->x *= rhs.x;
		this->y *= rhs.y;
		this->z *= rhs.z;
		return *this;
	}
	v3d_generic &operator/=(const v3d_generic &rhs)
	{
		this->x /= rhs.x;
		this->y /= rhs.y;
		this->z /= rhs.z;
		return *this;
	}
	v3d_generic operator+() const { return {+x, +y, +z}; }
	v3d_generic operator-() const { return {-x, -y, -z}; }
	bool operator==(const v3d_generic &rhs) const { return (this->x == rhs.x && this->y == rhs.y && this->z == rhs.z); }
	bool operator!=(const v3d_generic &rhs) const { return (this->x != rhs.x || this->y != rhs.y || this->z != rhs.z); }
	const std::string str() const { return std::string("(") + std::to_string(this->x) + "," + std::to_string(this->y) + "," + std::to_string(this->z) + ")"; }
	friend std::ostream &operator<<(std::ostream &os, const v3d_generic &rhs)
	{
		os << rhs.str();
		return os;
	}
	operator v3d_generic<int32_t>() const { return {static_cast<int32_t>(this->x), static_cast<int32_t>(this->y), static_cast<int32_t>(this->z)}; }
	operator v3d_generic<float>() const { return {static_cast<float>(this->x), static_cast<float>(this->y), static_cast<float>(this->z)}; }
	operator v3d_generic<double>() const { return {static_cast<double>(this->x), static_cast<double>(this->y), static_cast<double>(this->z)}; }

	//custom
	double sum() const { return this->x + this->y + this->z; }
	v3d_generic rotate(v3d_generic &a, double angle)
	{
		double c = std::cos(angle), s = std::sin(angle), cc = 1 - c;
		double xy = a.x * a.y, xz = a.x * a.z, yz = a.y * a.z;
		return {
			(c + cc * a.x * a.x) * x + (cc * xy - s * a.z) * y + (cc * xz + s * a.y) * z,
			(cc * xy + s * a.z) * x + (c + cc * a.y * a.y) * y + (cc * yz - s * a.x) * z,
			(cc * xz - s * a.y) * x + (cc * yz + s * a.x) * y + (c + cc * a.z * a.z) * z};
		//v3d_generic<T> k = a.norm();
		//double c = std::cos(angle), s = std::sin(angle);
		//return (*this * c) + k.cross(*this) * s + k * k.dot(*this) * (1 - c);
	}
	T max_value() const { return std::max(x, std::max(y, z)); }
	T min_value() const { return std::min(x, std::min(y, z)); }
};

typedef v3d_generic<int32_t> vi3d;
typedef v3d_generic<uint32_t> vu3d;
typedef v3d_generic<float> vf3d;
typedef v3d_generic<double> vd3d;

class GameObject;

struct svo_space
{
	std::vector<GameObject> items;
	std::vector<svo_space> childs;
	svo_space *parent = nullptr;
	uint32_t depth = 1;
	uint8_t used_childs_mask = 0x00;

	void set_used_child(svo_space *child)
	{
		if (child == &childs[0])
		{
			used_childs_mask |= 0x01;
			return;
		}
		if (child == &childs[1])
		{
			used_childs_mask |= 0x02;
			return;
		}
		if (child == &childs[2])
		{
			used_childs_mask |= 0x04;
			return;
		}
		if (child == &childs[3])
		{
			used_childs_mask |= 0x08;
			return;
		}
		if (child == &childs[4])
		{
			used_childs_mask |= 0x10;
			return;
		}
		if (child == &childs[5])
		{
			used_childs_mask |= 0x20;
			return;
		}
		if (child == &childs[6])
		{
			used_childs_mask |= 0x40;
			return;
		}
		if (child == &childs[7])
		{
			used_childs_mask |= 0x80;
			return;
		}
	}
	bool create_childs()
	{
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

struct material
{
	olc::Pixel color;
};

uint32_t nodes_count = 0;
uint32_t parent_count = 0;

struct svo_model
{
	const static uint32_t no_childs = 0x00ffffff;
	uint32_t material_id : 8, childs_id : 24;

	friend GameObject;

	svo_model() {
		nodes_count++;
		material_id = 0;
		childs_id = no_childs;
	}
	svo_model(const svo_model& model) {
		material_id = model.material_id;
		childs_id = model.childs_id;
	}
	svo_model(svo_model&& model) noexcept {
		material_id = model.material_id;
		childs_id = model.childs_id;
		model.childs_id = no_childs;
	}
	svo_model& operator=(const svo_model& model) {
		material_id = model.material_id;
		childs_id = model.childs_id;
		return *this;
	}
	svo_model& operator=(svo_model&& model) noexcept {
		material_id = model.material_id;
		childs_id = model.childs_id;
		model.childs_id = no_childs;
		return *this;
	}
	~svo_model() {}

	inline bool is_terminal() {
		return childs_id == no_childs;
	}


private:
	//bool create_childs(GameObject& go) {
	//	if (childs_id != no_childs)
	//		go.destroyChildNodes(childs_id);

	//	childs_id = go.createChildNodes();

	//	if (childs_id == no_childs)
	//		return false;

	//	parent_count++;

	//	return true;
	//}

	//void generator_random(GameObject& go) {
	//	if (max_depth == 1) {
	//		material_id = rand() % 255;
	//		return;
	//	}

	//	create_childs();
	//	uint8_t mask = rand();
	//	for (uint8_t i = 0; i < 8; i++) {
	//		if (mask & (1 << i)) {
	//			childs[i].generator_random();
	//		}
	//		else {
	//			childs[i].material_id = 255;
	//		}
	//	}
	//}

	//void generator_total_split(GameObject& go) {
	//	if (max_depth == 1) {
	//		material_id = rand() % 255;
	//		return;
	//	}

	//	create_childs();
	//	for (int i = 0; i < 8; i++) {
	//		childs[i].generator_total_split();
	//	}
	//}

	//void generator_random_2(GameObject& go) {
	//	if (max_depth == 1) {
	//		material_id = rand() % 255;
	//		return;
	//	}

	//	if (!(rand() & 3)) {
	//		material_id = rand() % 255;
	//		return;
	//	}

	//	create_childs();
	//	uint8_t mask = rand();
	//	for (uint8_t i = 0; i < 8; i++) {
	//		if (mask & (1 << i)) {
	//			childs[i].generator_random_2();
	//		}
	//		else {
	//			childs[i].material_id = 255;
	//		}
	//	}
	//}

	//void generator_sphere_corner(GameObject& go) {
	//	if (max_depth == 1) {
	//		material_id = rand() % 255;
	//		return;
	//	}

	//	create_childs();
	//	childs[0]._generator_sphere(0.0, 0.0, 0.0, 0.5, 0.5, 0.5);
	//	childs[1]._generator_sphere(0.5, 0.0, 0.0, 1.0, 0.5, 0.5);
	//	childs[2]._generator_sphere(0.0, 0.5, 0.0, 0.5, 1.0, 0.5);
	//	childs[3]._generator_sphere(0.5, 0.5, 0.0, 1.0, 1.0, 0.5);
	//	childs[4]._generator_sphere(0.0, 0.0, 0.5, 0.5, 0.5, 1.0);
	//	childs[5]._generator_sphere(0.5, 0.0, 0.5, 1.0, 0.5, 1.0);
	//	childs[6]._generator_sphere(0.0, 0.5, 0.5, 0.5, 1.0, 1.0);
	//	childs[7]._generator_sphere(0.5, 0.5, 0.5, 1.0, 1.0, 1.0);
	//}
};

// WIP
struct svo_child_wrapper
{
	svo_model childs[8];
};

class GameObject
{
private:
	uint8_t max_depth = 10;
	std::vector<svo_model*> child_nodes;
	std::vector<uint32_t> unused_nodes;


public:
	const static uint32_t MAX_NODES = 0x00ffffff;
	svo_model top_node;
	vd3d size;
	vd3d pos;
	material materials[256];


public:
	GameObject() {
		for (int i = 0; i < 255; i++)
			materials[i].color = olc::Pixel(rand(), rand(), rand(), 255);
		materials[255].color = olc::BLANK;
	}
	GameObject(const GameObject& go)
		: top_node(go.top_node)
	{
		for (auto& c : go.child_nodes)
			child_nodes.push_back(c);
		for (int i = 0; i < 256; i++)
			materials[i] = go.materials[i];
	}


public:
	svo_model* getChilds(uint32_t child_id) {
		return child_nodes[child_id];
	}


	bool generateModelSphere(uint8_t depth) {
		if (depth == 0)
			return false;

		clear();

		depth = std::min(depth, max_depth);

		if (!_generateModelSphere(
			depth, top_node,
			-1.0, -1.0, -1.0,
			1.0, 1.0, 1.0
		)) {
			clear();
			return false;
		}
		return true;
	}
	bool generateModelSlope(uint8_t depth) {
		if (depth == 0)
			return false;

		clear();

		depth = std::min(depth, max_depth);

		if (!_generateModelSlope(depth, top_node)) {
			clear();
			return false;
		}
		return true;
	}
	bool generateModelTotalSplit(uint8_t depth) {
		if (depth == 0)
			return false;

		clear();

		depth = std::min(depth, max_depth);

		if (!_generateModelTotalSplit(depth, top_node)) {
			clear();
			return false;
		}
		return true;
	}
	//bool generateModelSerpinsky(uint8_t depth) {
	//	if (depth == 0)
	//		return false;

	//	clear();

	//	depth = std::min(depth, max_depth);

	//	if (!_generateModelSerpinsky()) {
	//		clear();
	//		return false;
	//	}
	//	return true;
	//}


	void clear() {
		for (auto c : child_nodes)
			delete[] c;
		unused_nodes.clear();
		child_nodes.clear();
	}


private:
	uint32_t createChildNodes() {
		if (!unused_nodes.empty()) {
			uint32_t nodes = unused_nodes.back();
			unused_nodes.pop_back();
			return nodes;
		}

		if (child_nodes.size() == MAX_NODES)
			return MAX_NODES;
		parent_count++;
		child_nodes.emplace_back(new svo_model[8]);
		return child_nodes.size() - 1;
	}

	bool _generateModelSphere(
		uint8_t depth, svo_model& model,
		double x0, double y0, double z0,
		double x1, double y1, double z1
	) {
		double xm = (x0 + x1) * 0.5;
		double ym = (y0 + y1) * 0.5;
		double zm = (z0 + z1) * 0.5;

		double xs = (x1 - x0);
		double ys = (y1 - y0);
		double zs = (z1 - z0);

		double center = std::sqrt(xm * xm + ym * ym + zm * zm);
		double diag = std::sqrt(xs * xs + ys * ys + zs * zs) * 0.5;

		if (depth == 1) {
			if (center < 1.0)
				model.material_id = rand() % 255;
			else
				model.material_id = 255;
			return true;
		}

		// full inside
		if (center + diag <= 1.0) {
			model.material_id = rand() % 255;
			return true;
		}

		// full outside
		if (center - diag >= 1.0) {
			model.material_id = 255;
			return true;
		}

		model.childs_id = createChildNodes();
		if (model.childs_id == MAX_NODES)
			return false;

		svo_model* childs = child_nodes[model.childs_id];
		if (
			!_generateModelSphere(depth - 1, childs[0], x0, y0, z0, xm, ym, zm) ||
			!_generateModelSphere(depth - 1, childs[1], xm, y0, z0, x1, ym, zm) ||
			!_generateModelSphere(depth - 1, childs[2], x0, ym, z0, xm, y1, zm) ||
			!_generateModelSphere(depth - 1, childs[3], xm, ym, z0, x1, y1, zm) ||
			!_generateModelSphere(depth - 1, childs[4], x0, y0, zm, xm, ym, z1) ||
			!_generateModelSphere(depth - 1, childs[5], xm, y0, zm, x1, ym, z1) ||
			!_generateModelSphere(depth - 1, childs[6], x0, ym, zm, xm, y1, z1) ||
			!_generateModelSphere(depth - 1, childs[7], xm, ym, zm, x1, y1, z1)
		) {
			return false;
		}

		model.material_id = childs[0].material_id;

		return true;
	}
	bool _generateModelSlope(uint8_t depth, svo_model& model) {
		if (depth == 1) {
			model.material_id = rand() % 255;
			return true;
		}

		model.childs_id = createChildNodes();
		if (model.childs_id == MAX_NODES)
			return false;

		svo_model* childs = child_nodes[model.childs_id];
		_generateModelSlope(depth - 1, childs[0]);
		_generateModelSlope(depth - 1, childs[3]);
		_generateModelSlope(depth - 1, childs[6]);
		childs[2].material_id = 255;
		childs[1].material_id = rand() % 255;
		childs[4].material_id = rand() % 255;
		childs[5].material_id = rand() % 255;
		childs[7].material_id = rand() % 255;

		model.material_id = childs[0].material_id;

		return true;
	}
	bool _generateModelTotalSplit(uint8_t depth, svo_model& model) {
		if (depth == 1) {
			model.material_id = rand() % 255;
			return true;
		}

		model.childs_id = createChildNodes();
		if (model.childs_id == MAX_NODES)
			return false;

		svo_model* childs = child_nodes[model.childs_id];
		if (
			!_generateModelTotalSplit(depth - 1, childs[0]) ||
			!_generateModelTotalSplit(depth - 1, childs[1]) ||
			!_generateModelTotalSplit(depth - 1, childs[2]) ||
			!_generateModelTotalSplit(depth - 1, childs[3]) ||
			!_generateModelTotalSplit(depth - 1, childs[4]) ||
			!_generateModelTotalSplit(depth - 1, childs[5]) ||
			!_generateModelTotalSplit(depth - 1, childs[6]) ||
			!_generateModelTotalSplit(depth - 1, childs[7])
		) {
			return false;
		}

		model.material_id = childs[0].material_id;

		return true;
	}
};

class Space3D
{
public:
	svo_space svo;
	vu3d size;
	vi3d pos;
	vi3d oppositePos;

public:
	Space3D() {}
	Space3D(vi3d &&coords, vu3d &&sizes)
	{
		size = sizes;
		pos = coords;
		oppositePos = pos + size;

		uint32_t maxSide = std::max(size.x, std::max(size.y, size.z));
		uint32_t depth = 1;
		uint32_t boxSide = 1;
		while (boxSide < maxSide)
		{
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

class Camera
{

public:
	vd3d pos;
	double hAngle, vAngle;
	double move_speed = 0.1;
	double turn_speed = 0.02;
	const double fovH = M_PI / 2; // 90 degree
	const double fovV = M_PI / 3; // 60 degree

public:
	Camera()
		: pos(0.0, 0.0, 0.0), hAngle(0.0), vAngle(0.0) {}
	Camera(double _x, double _y, double _z)
		: pos(_x, _y, _z), hAngle(0.0), vAngle(0.0) {}
	Camera(double _x, double _y, double _z, double _hAngle, double _vAngle)
		: pos(_x, _y, _z), hAngle(_hAngle), vAngle(_vAngle) {}

public:
	void turnRight()
	{
		hAngle -= turn_speed;
		if (hAngle < -M_PI)
			hAngle += M_PI + M_PI;
	}
	void turnLeft()
	{
		hAngle += turn_speed;
		if (hAngle > M_PI)
			hAngle -= M_PI + M_PI;
	}
	void turnUp()
	{
		if (vAngle <= M_PI / 2.0)
			vAngle += turn_speed;
	}
	void turnDown()
	{
		if (vAngle >= -M_PI / 2.0)
			vAngle -= turn_speed;
	}

	void moveW()
	{
		pos.x += move_speed * std::cos(hAngle) * std::cos(vAngle);
		pos.y += move_speed * std::sin(vAngle);
		pos.z += move_speed * -std::sin(hAngle) * std::cos(vAngle);
	}
	void moveA()
	{
		double _hAngle = hAngle + M_PI / 2.0;
		pos.x += move_speed * std::cos(_hAngle) * std::cos(vAngle);
		pos.y += move_speed * std::sin(vAngle);
		pos.z += move_speed * -std::sin(_hAngle) * std::cos(vAngle);
	}
	void moveS()
	{
		pos.x -= move_speed * std::cos(hAngle) * std::cos(vAngle);
		pos.y -= move_speed * std::sin(vAngle);
		pos.z -= move_speed * -std::sin(hAngle) * std::cos(vAngle);
	}
	void moveD()
	{
		double _hAngle = hAngle - M_PI / 2.0;
		pos.x += move_speed * std::cos(_hAngle) * std::cos(vAngle);
		pos.y += move_speed * std::sin(vAngle);
		pos.z += move_speed * -std::sin(_hAngle) * std::cos(vAngle);
	}
	void moveUp()
	{
		pos.x += move_speed * -std::sin(vAngle) * -std::cos(hAngle);
		pos.y += move_speed * std::cos(vAngle);
		pos.z += move_speed * -std::sin(vAngle) * -std::sin(hAngle);
	}
	void moveDown()
	{
		pos.x -= move_speed * -std::sin(vAngle) * -std::cos(hAngle);
		pos.y -= move_speed * std::cos(vAngle);
		pos.z -= move_speed * -std::sin(vAngle) * -std::sin(hAngle);
	}

	vd3d norm_direction() const
	{
		return {
			std::cos(hAngle) * std::cos(vAngle),
			std::sin(vAngle),
			-std::sin(hAngle) * std::cos(vAngle),
		};
	}
	vd3d perp_vert() const
	{
		double _vAngle = vAngle + M_PI / 2.0;
		return {
			std::cos(hAngle) * std::cos(_vAngle),
			std::sin(_vAngle),
			-std::sin(hAngle) * std::cos(_vAngle),
		};
	}
	vd3d perp_hor() const
	{
		double _hAngle = hAngle + M_PI / 2.0;
		return {
			-std::cos(_hAngle) * std::cos(vAngle),
			std::sin(vAngle),
			std::sin(_hAngle) * std::cos(vAngle),
		};
	}
};


struct worker {
	bool alive = true;
	bool started = false;
	std::mutex mut;
	std::thread thread;

	vd3d corner;
	vd3d deltaX, deltaY;
	std::function<void(vd3d, vd3d, vd3d)> _job;


	worker() {}
	worker(const worker& w) = delete;
	worker(worker&& w) noexcept : thread(std::move(w.thread)) {}


	void job(std::condition_variable& work_cv, std::condition_variable& complete_cv) {
		std::unique_lock<std::mutex> ul(mut);
		started = true;

		work_cv.wait(ul);
		while (alive) {
			_job(corner, deltaX, deltaY);
			complete_cv.notify_one();
			work_cv.wait(ul);
		}
	}

	void set_params(vd3d _luCorner, vd3d _deltaX, vd3d _deltaY) {
		corner = _luCorner;
		deltaX = _deltaX;
		deltaY = _deltaY;
	}
};


enum class axis : uint8_t {
	None = 0,
	X = 1,
	Y = 2,
	Z = 3,
};

std::ostream& operator<<(std::ostream& os, const axis& axis) {
	switch (axis)
	{
	case axis::X:
		os << "X";
		break;
	case axis::Y:
		os << "Y";
		break;
	case axis::Z:
		os << "Z";
		break;
	case axis::None:
	default:
		os << "None";
		break;
	}

	return os;
}

struct ray_render_info {
	vd3d ray_source = {};
	vd3d ray_direction = {};
	svo_model* model = nullptr;
	uint8_t depth = 0;
	uint8_t max_depth = 0;
	uint16_t func_calls = 0;
	double axis_param = 0.0;
	double param = 0.0;
	axis axis = axis::None;

	friend std::ostream& operator<<(std::ostream& os, const ray_render_info& rri);
};

std::ostream& operator<<(std::ostream& os, const ray_render_info& rri) {
	os << "ray_source: " << rri.ray_source << '\n'
		<< "ray_direction: " << rri.ray_direction << '\n'
		<< "depth: " << (uint32_t)rri.depth << '\n'
		<< "max_depth: " << (uint32_t)rri.max_depth << '\n'
		<< "func_calls: " << (uint32_t)rri.func_calls << '\n'
		<< "axis_param: " << rri.axis_param << '\n'
		<< "param: " << rri.param << '\n'
		<< "axis: " << rri.axis << '\n';
	return os;
}


struct config {
	double min_view_size = 0.00001;
};


class SVOGame;

typedef void (SVOGame::* RenderFunc)(uint32_t, uint32_t, GameObject&, ray_render_info&);


class SVOGame : public olc::PixelGameEngine
{
public:
	SVOGame()
	{
		sAppName = "SVO";
	}


protected:
	Camera player_view = Camera(-5.0, 0.0, 0.0);
	Space3D space;
	GameObject testObject;
	config config;

	uint32_t stop_x, stop_y;

	std::mutex mut;
	std::condition_variable work_cv;
	std::condition_variable complete_cv;
	std::vector<worker> workers;
	std::atomic_uint8_t work_in_progress;

	std::vector<RenderFunc> render_funcs;
	std::vector<RenderFunc>::iterator render_func_it;


protected:
	bool OnUserCreate() override
	{
		svo_model model;
		srand(time(nullptr));

		stop_x = ScreenWidth();
		stop_y = ScreenHeight();
		
		uint8_t depth = 7;
		double size = 16.0;
		testObject.pos = { 0.0, -size / 2.0, -size / 2.0 };
		testObject.size = { size, size, size };
		if (!testObject.generateModelSphere(depth))
			std::cout << "Error on creating model" << std::endl;
		

		render_funcs.push_back(&SVOGame::RenderDefault);
		render_funcs.push_back(&SVOGame::RenderLength);
		render_funcs.push_back(&SVOGame::RenderDepth);
		render_funcs.push_back(&SVOGame::RenderMaxDepth);
		render_funcs.push_back(&SVOGame::RenderFuncCalls);

		render_func_it = render_funcs.begin();

		setWorkerThreads(1);

		std::cout << "nodes_count: " << nodes_count << std::endl;
		std::cout << "parent_count: " << parent_count << std::endl;

		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		std::unique_lock<std::mutex> ul(mut);

		// Controls
		if (GetKey(olc::ESCAPE).bPressed)
			return false;

		if (GetKey(olc::CTRL).bHeld) {
			if (GetKey(olc::EQUALS).bPressed) {
				uint8_t count = workers.size();
				if (count < 32)
					setWorkerThreads(count + 1);
			}
			if (GetKey(olc::MINUS).bPressed) {
				uint8_t count = workers.size();
				if (count > 1)
					setWorkerThreads(count - 1);
			}
		}
		else {
			if (GetKey(olc::MINUS).bPressed)
				player_view.move_speed /= 2.0;
			if (GetKey(olc::EQUALS).bPressed)
				player_view.move_speed *= 2.0;
		}

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

		if (GetKey(olc::R).bPressed) {
			render_func_it++;
			if (render_func_it == render_funcs.end())
				render_func_it = render_funcs.begin();
		}

		if (GetKey(olc::T).bPressed)
			config.min_view_size -= 0.000001;
		if (GetKey(olc::G).bPressed)
			config.min_view_size += 0.000001;

		if (GetMouse(0).bPressed) {
			stop_x = GetMouseX();
			stop_y = GetMouseY();
		}

		// Render
		Clear(olc::BLACK);

		const double multiplier = 1000.0;

		vd3d vert = player_view.perp_vert() * std::tan(player_view.fovV / 2.0) / multiplier;
		vd3d hor = player_view.perp_hor() * std::tan(player_view.fovH / 2.0) / multiplier;

		vd3d norm_dir = player_view.norm_direction() / multiplier;

		vd3d luCorner = player_view.pos + norm_dir + vert - hor;
		vd3d ruCorner = player_view.pos + norm_dir + vert + hor;
		vd3d lbCorner = player_view.pos + norm_dir - vert - hor;
		vd3d rbCorner = player_view.pos + norm_dir - vert + hor;

		vd3d deltaY = (lbCorner - luCorner) / ScreenHeight();
		vd3d deltaX = (ruCorner - luCorner) / ScreenWidth();

		for (int i = 0; i < workers.size(); i++) {
			workers[i].set_params(luCorner, deltaX, deltaY);
		}

		work_in_progress = workers.size();
		work_cv.notify_all();
		while (work_in_progress > 0) { std::this_thread::yield(); }
		//complete_cv.wait(ul, [&] { return work_in_progress == 0; });

		olc::vi2d mousePos = GetMousePos();

		std::stringstream coords;
		coords << "pos: " << player_view.pos.str() << '\n';
		coords << "dir: " << player_view.norm_direction().str() << '\n';
		coords << "mouse dir: " << (luCorner + deltaY * mousePos.y + deltaX * mousePos.x - player_view.pos).str() << '\n';
		coords << "move_speed: " << player_view.move_speed << '\n';
		coords << "worker threads: " << workers.size() << '\n';
		coords << "config.min_view_size: " << config.min_view_size << '\n';
		DrawString(0, 0, coords.str());

		return true;
	}

	bool OnUserDestroy() override {
		setWorkerThreads(0);
		testObject.clear();
		return true;
	}

	
protected:
	ray_render_info rayParameterV2(GameObject& go, vd3d ray_source, vd3d ray_dir) {
		unsigned char a = 0;
		ray_render_info rri;

		rri.ray_source = ray_source;
		rri.ray_direction = ray_dir;

		ray_source -= go.pos;

		if (ray_dir.x < 0.0)
		{
			ray_source.x = go.size.x - ray_source.x;
			ray_dir.x = -ray_dir.x;
			a |= 1;
		}
		if (ray_dir.y < 0.0)
		{
			ray_source.y = go.size.y - ray_source.y;
			ray_dir.y = -ray_dir.y;
			a |= 2;
		}
		if (ray_dir.z < 0.0)
		{
			ray_source.z = go.size.z - ray_source.z;
			ray_dir.z = -ray_dir.z;
			a |= 4;
		}

		if (ray_dir.x == 0.0) {
			if (ray_source.x < 0.0 || ray_source.x >= go.size.x)
				return rri;
			ray_dir.x = 0.0000000001;
		}
		if (ray_dir.y == 0.0) {
			if (ray_source.y < 0.0 || ray_source.y >= go.size.y)
				return rri;
			ray_dir.y = 0.0000000001;
		}
		if (ray_dir.z == 0.0) {
			if (ray_source.z < 0.0 || ray_source.z >= go.size.z)
				return rri;
			ray_dir.z = 0.0000000001;
		}

		double divX = 1.0 / ray_dir.x;
		double divY = 1.0 / ray_dir.y;
		double divZ = 1.0 / ray_dir.z;

		double tx0 = -ray_source.x * divX;
		double tx1 = (go.size.x - ray_source.x) * divX;
		double ty0 = -ray_source.y * divY;
		double ty1 = (go.size.y - ray_source.y) * divY;
		double tz0 = -ray_source.z * divZ;
		double tz1 = (go.size.z - ray_source.z) * divZ;

		double txm = (tx0 + tx1) * 0.5;
		double tym = (ty0 + ty1) * 0.5;
		double tzm = (tz0 + tz1) * 0.5;

		if (std::max(tx0, std::max(ty0, tz0)) >= std::min(tx1, std::min(ty1, tz1)))
			return rri;

		if (tx1 < 0.0 || ty1 < 0.0 || tz1 < 0.0)
			return rri;

		struct stack_struct {
			vd3d t0;
			vd3d tm;
			vd3d t1;
			vd3d node_size;
			svo_model* node = nullptr;
			uint8_t child_id = 0;
		};
		stack_struct stack[11];
		int32_t stack_idx = 0;
		uint8_t depth = 1;
		stack[0] = {
			{ tx0, ty0, tz0 },
			{ txm, tym, tzm },
			{ tx1, ty1, tz1 },
			go.size,
			&go.top_node,
			firstNodeXYZ(tx0, ty0, tz0, txm, tym, tzm),
		};

		auto nextNode = [this](uint8_t current, vd3d& tm, vd3d& t1) {
			switch (current) {
			case 0:
				return nextNode3D(tm.x, 1, tm.y, 2, tm.z, 4);
			case 1:
				return nextNode3D(t1.x, 8, tm.y, 3, tm.z, 5);
			case 2:
				return nextNode3D(tm.x, 3, t1.y, 8, tm.z, 6);
			case 3:
				return nextNode3D(t1.x, 8, t1.y, 8, tm.z, 7);
			case 4:
				return nextNode3D(tm.x, 5, tm.y, 6, t1.z, 8);
			case 5:
				return nextNode3D(t1.x, 8, tm.y, 7, t1.z, 8);
			case 6:
				return nextNode3D(tm.x, 7, t1.y, 8, t1.z, 8);
			case 7:
				return nextNode3D(t1.x, 8, t1.y, 8, t1.z, 8);
			default:
				return uint8_t(8);
			}
		};

		while (stack_idx >= 0) {
			rri.max_depth = std::max(rri.max_depth, depth);
			rri.func_calls++;

			stack_struct& current_stack = stack[stack_idx];
			vd3d& t0 = current_stack.t0;
			vd3d& tm = current_stack.tm;
			vd3d& t1 = current_stack.t1;
			vd3d& node_size = current_stack.node_size;
			svo_model* node = current_stack.node;
			uint8_t& child_id = current_stack.child_id;

			if (child_id == 8) {
				stack_idx--;
				depth--;
				continue;
			}

			svo_model* child = go.getChilds(node->childs_id) + (child_id ^ a);
			vd3d child_t0 = {
				child_id & 0x01 ? tm.x : t0.x,
				child_id & 0x02 ? tm.y : t0.y,
				child_id & 0x04 ? tm.z : t0.z,
			};
			vd3d child_t1 = {
				child_id & 0x01 ? t1.x : tm.x,
				child_id & 0x02 ? t1.y : tm.y,
				child_id & 0x04 ? t1.z : tm.z,
			};
			vd3d child_tm = (child_t0 + child_t1) * 0.5;

			if (child_t1.x < 0.0 || child_t1.y < 0.0 || child_t1.z < 0.0) {
				child_id = nextNode(child_id, tm, t1);
				if (child_id == 8) {
					stack_idx--;
					depth--;
				}
				continue;
			}
			double param = child_t0.max_value();
			if (
				child->is_terminal() ||
				std::abs(node_size.sum() / param) < config.min_view_size
			) {
				if (go.materials[child->material_id].color.a != 255) {
					child_id = nextNode(child_id, tm, t1);
					if (child_id == 8) {
						stack_idx--;
						depth--;
					}
					continue;
				}
				rri.model = child;
				rri.depth = depth;
				rri.param = param;
				// max
				if (t0.x > t0.y) {
					if (t0.x > t0.z) {
						rri.axis_param = rri.ray_direction.x;
						rri.axis = axis::X;
					}
					else {
						rri.axis_param = rri.ray_direction.z;
						rri.axis = axis::Z;
					}
				}
				else if (t0.y > t0.z) {
					rri.axis_param = rri.ray_direction.y;
					rri.axis = axis::Y;
				}
				else {
					rri.axis_param = rri.ray_direction.z;
					rri.axis = axis::Z;
				}
				break;
			}

			child_id = nextNode(child_id, tm, t1);
			uint8_t next_child_id = firstNodeXYZ(
				child_t0.x, child_t0.y, child_t0.z,
				child_tm.x, child_tm.y, child_tm.z
			);
			stack_idx++;
			depth++;
			stack[stack_idx] = {
				child_t0,
				child_tm,
				child_t1,
				node_size * 0.5,
				child,
				next_child_id,
			};
		}
		return rri;
	}

	ray_render_info rayParameter(GameObject& go, vd3d ray_source, vd3d ray_dir)
	{
		unsigned char a = 0;
		ray_render_info rri;

		rri.ray_source = ray_source;
		rri.ray_direction = ray_dir;

		ray_source -= go.pos;

		if (ray_dir.x < 0.0)
		{
			ray_source.x = go.size.x - ray_source.x;
			ray_dir.x = -ray_dir.x;
			a |= 1;
		}
		if (ray_dir.y < 0.0)
		{
			ray_source.y = go.size.y - ray_source.y;
			ray_dir.y = -ray_dir.y;
			a |= 2;
		}
		if (ray_dir.z < 0.0)
		{
			ray_source.z = go.size.z - ray_source.z;
			ray_dir.z = -ray_dir.z;
			a |= 4;
		}

		if (ray_dir.x == 0.0) {
			if (ray_source.x < 0.0 || ray_source.x >= go.size.x)
				return rri;

			if (ray_dir.y == 0.0) {
				if (ray_source.y < 0.0 || ray_source.y >= go.size.y)
					return rri;
				double divZ = 1.0 / ray_dir.z;
				double tz0 = -ray_source.z * divZ;
				double tz1 = (go.size.z - ray_source.z) * divZ;
				procSubtreeZ(tz0, tz1, go, &go.top_node, 1, a, go.size, ray_source, rri);
				return rri;
			}
			else if (ray_dir.z == 0.0) {
				if (ray_source.z < 0.0 || ray_source.z >= go.size.z)
					return rri;
				double divY = 1.0 / ray_dir.y;
				double ty0 = -ray_source.y * divY;
				double ty1 = (go.size.y - ray_source.y) * divY;
				procSubtreeY(ty0, ty1, go, &go.top_node, 1, a, go.size, ray_source, rri);
				return rri;
			}
			double divY = 1.0 / ray_dir.y;
			double divZ = 1.0 / ray_dir.z;
			double ty0 = -ray_source.y * divY;
			double ty1 = (go.size.y - ray_source.y) * divY;
			double tz0 = -ray_source.z * divZ;
			double tz1 = (go.size.z - ray_source.z) * divZ;
			if (std::max(ty0, tz0) < std::min(ty1, tz1))
				procSubtreeYZ(ty0, tz0, ty1, tz1, go, &go.top_node, 1, a, go.size, ray_source, rri);
			return rri;
		}
		if (ray_dir.y == 0.0) {
			if (ray_source.y < 0.0 || ray_source.y >= go.size.y)
				return rri;
			if (ray_dir.z == 0.0) {
				if (ray_source.z < 0.0 || ray_source.z >= go.size.z)
					return rri;
				double divX = 1.0 / ray_dir.x;
				double tx0 = -ray_source.x * divX;
				double tx1 = (go.size.x - ray_source.x) * divX;
				procSubtreeX(tx0, tx1, go, &go.top_node, 1, a, go.size, ray_source, rri);
				return rri;
			}
			double divX = 1.0 / ray_dir.x;
			double divZ = 1.0 / ray_dir.z;
			double tx0 = -ray_source.x * divX;
			double tx1 = (go.size.x - ray_source.x) * divX;
			double tz0 = -ray_source.z * divZ;
			double tz1 = (go.size.z - ray_source.z) * divZ;
			if (std::max(tx0, tz0) < std::min(tx1, tz1))
				procSubtreeXZ(tx0, tz0, tx1, tz1, go, &go.top_node, 1, a, go.size, ray_source, rri);
			return rri;
		}
		if (ray_dir.z == 0.0) {
			if (ray_source.z < 0.0 || ray_source.z >= go.size.z)
				return rri;
			double divX = 1.0 / ray_dir.x;
			double divY = 1.0 / ray_dir.y;
			double tx0 = -ray_source.x * divX;
			double tx1 = (go.size.x - ray_source.x) * divX;
			double ty0 = -ray_source.y * divY;
			double ty1 = (go.size.y - ray_source.y) * divY;
			if (std::max(tx0, ty0) < std::min(tx1, ty1))
				procSubtreeXY(tx0, ty0, tx1, ty1, go, &go.top_node, 1, a, go.size, ray_source, rri);
			return rri;
		}

		double divX = 1.0 / ray_dir.x;
		double divY = 1.0 / ray_dir.y;
		double divZ = 1.0 / ray_dir.z;

		double tx0 = -ray_source.x * divX;
		double tx1 = (go.size.x - ray_source.x) * divX;
		double ty0 = -ray_source.y * divY;
		double ty1 = (go.size.y - ray_source.y) * divY;
		double tz0 = -ray_source.z * divZ;
		double tz1 = (go.size.z - ray_source.z) * divZ;

		if (std::max(tx0, std::max(ty0, tz0)) < std::min(tx1, std::min(ty1, tz1)))
			procSubtreeXYZ(tx0, ty0, tz0, tx1, ty1, tz1, go, &go.top_node, 1, a, go.size, rri);
		return rri;
	}

	void procSubtreeX(double tx0, double tx1, GameObject& go, svo_model* node, uint8_t depth, uint8_t a, vd3d node_1, vd3d ray_source, ray_render_info& ray_info) {
		ray_info.max_depth = std::max(ray_info.max_depth, depth);
		ray_info.func_calls++;

		if (tx1 < 0.0)
			return;
		double param = tx0;
		if (
			node->is_terminal() ||
			std::abs(node_1.sum() / param) < config.min_view_size
		) {
			if (go.materials[node->material_id].color.a != 255)
				return;
			ray_info.model = node;
			ray_info.depth = depth;
			ray_info.param = param;
			ray_info.axis_param = std::abs(ray_info.ray_direction.x);
			ray_info.axis = axis::X;
			return;
		}

		double txm = 0.5 * (tx0 + tx1);

		node_1 *= 0.5;

		uint8_t curr_node = firstNodeX(ray_source, node_1);

		do {
			switch (curr_node) {
			case 0:
				procSubtreeX(tx0, txm, go, go.getChilds(node->childs_id) + (a), depth + 1, a, node_1, ray_source, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode1D(1);
				break;
			case 1:
				procSubtreeX(txm, tx1, go, go.getChilds(node->childs_id) + (a ^ 1), depth + 1, a, node_1, ray_source - vd3d{ node_1.x, 0.0, 0.0 }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode1D(8);
				break;
			case 2:
				procSubtreeX(tx0, txm, go, go.getChilds(node->childs_id) + (a ^ 2), depth + 1, a, node_1, ray_source - vd3d{ 0.0, node_1.y, 0.0 }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode1D(3);
				break;
			case 3:
				procSubtreeX(txm, tx1, go, go.getChilds(node->childs_id) + (a ^ 3), depth + 1, a, node_1, ray_source - vd3d{ node_1.x, node_1.y, 0.0 }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode1D(8);
				break;
			case 4:
				procSubtreeX(tx0, txm, go, go.getChilds(node->childs_id) + (a ^ 4), depth + 1, a, node_1, ray_source - vd3d{ 0.0, 0.0, node_1.z }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode1D(5);
				break;
			case 5:
				procSubtreeX(txm, tx1, go, go.getChilds(node->childs_id) + (a ^ 5), depth + 1, a, node_1, ray_source - vd3d{ node_1.x, 0.0, node_1.z }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode1D(8);
				break;
			case 6:
				procSubtreeX(tx0, txm, go, go.getChilds(node->childs_id) + (a ^ 6), depth + 1, a, node_1, ray_source - vd3d{ 0.0, node_1.y, node_1.z }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode1D(7);
				break;
			case 7:
				procSubtreeX(txm, tx1, go, go.getChilds(node->childs_id) + (a ^ 7), depth + 1, a, node_1, ray_source - vd3d{ node_1.x, node_1.y, node_1.z }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode1D(8);
				break;
			}
		} while (curr_node < 8);
	}

	void procSubtreeY(double ty0, double ty1, GameObject& go, svo_model* node, uint8_t depth, uint8_t a, vd3d node_1, vd3d ray_source, ray_render_info& ray_info) {
		ray_info.max_depth = std::max(ray_info.max_depth, depth);
		ray_info.func_calls++;

		if (ty1 < 0.0)
			return;
		double param = ty0;
		if (
			node->is_terminal() ||
			std::abs(node_1.sum() / param) < config.min_view_size
		) {
			if (go.materials[node->material_id].color.a != 255)
				return;
			ray_info.model = node;
			ray_info.depth = depth;
			ray_info.param = param;
			ray_info.axis_param = std::abs(ray_info.ray_direction.y);
			ray_info.axis = axis::Y;
			return;
		}

		double tym = 0.5 * (ty0 + ty1);

		node_1 *= 0.5;

		uint8_t curr_node = firstNodeY(ray_source, node_1);

		do {
			switch (curr_node) {
			case 0:
				procSubtreeY(ty0, tym, go, go.getChilds(node->childs_id) + (a), depth + 1, a, node_1, ray_source, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode1D(2);
				break;
			case 1:
				procSubtreeY(ty0, tym, go, go.getChilds(node->childs_id) + (a ^ 1), depth + 1, a, node_1, ray_source - vd3d{ node_1.x, 0.0, 0.0 }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode1D(3);
				break;
			case 2:
				procSubtreeY(tym, ty1, go, go.getChilds(node->childs_id) + (a ^ 2), depth + 1, a, node_1, ray_source - vd3d{ 0.0, node_1.y, 0.0 }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode1D(8);
				break;
			case 3:
				procSubtreeY(tym, ty1, go, go.getChilds(node->childs_id) + (a ^ 3), depth + 1, a, node_1, ray_source - vd3d{ node_1.x, node_1.y, 0.0 }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode1D(8);
				break;
			case 4:
				procSubtreeY(ty0, tym, go, go.getChilds(node->childs_id) + (a ^ 4), depth + 1, a, node_1, ray_source - vd3d{ 0.0, 0.0, node_1.z }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode1D(6);
				break;
			case 5:
				procSubtreeY(ty0, tym, go, go.getChilds(node->childs_id) + (a ^ 5), depth + 1, a, node_1, ray_source - vd3d{ node_1.x, 0.0, node_1.z }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode1D(7);
				break;
			case 6:
				procSubtreeY(tym, ty1, go, go.getChilds(node->childs_id) + (a ^ 6), depth + 1, a, node_1, ray_source - vd3d{ 0.0, node_1.y, node_1.z }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode1D(8);
				break;
			case 7:
				procSubtreeY(tym, ty1, go, go.getChilds(node->childs_id) + (a ^ 7), depth + 1, a, node_1, ray_source - vd3d{ node_1.x, node_1.y, node_1.z }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode1D(8);
				break;
			}
		} while (curr_node < 8);
	}

	void procSubtreeZ(double tz0, double tz1, GameObject& go, svo_model* node, uint8_t depth, uint8_t a, vd3d node_1, vd3d ray_source, ray_render_info& ray_info) {
		ray_info.max_depth = std::max(ray_info.max_depth, depth);
		ray_info.func_calls++;

		if (tz1 < 0.0)
			return;
		double param = tz0;
		if (
			node->is_terminal() ||
			std::abs(node_1.sum() / param) < config.min_view_size
		) {
			if (go.materials[node->material_id].color.a != 255)
				return;
			ray_info.model = node;
			ray_info.depth = depth;
			ray_info.param = param;
			ray_info.axis_param = std::abs(ray_info.ray_direction.z);
			ray_info.axis = axis::Z;
			return;
		}

		double tzm = 0.5 * (tz0 + tz1);

		node_1 *= 0.5;

		uint8_t curr_node = firstNodeZ(ray_source, node_1);

		do {
			switch (curr_node) {
			case 0:
				procSubtreeZ(tz0, tzm, go, go.getChilds(node->childs_id) + (a), depth + 1, a, node_1, ray_source, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode1D(4);
				break;
			case 1:
				procSubtreeZ(tz0, tzm, go, go.getChilds(node->childs_id) + (a ^ 1), depth + 1, a, node_1, ray_source - vd3d{ node_1.x, 0.0, 0.0 }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode1D(5);
				break;
			case 2:
				procSubtreeZ(tz0, tzm, go, go.getChilds(node->childs_id) + (a ^ 2), depth + 1, a, node_1, ray_source - vd3d{ 0.0, node_1.y, 0.0 }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode1D(6);
				break;
			case 3:
				procSubtreeZ(tz0, tzm, go, go.getChilds(node->childs_id) + (a ^ 3), depth + 1, a, node_1, ray_source - vd3d{ node_1.x, node_1.y, 0.0 }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode1D(7);
				break;
			case 4:
				procSubtreeZ(tzm, tz1, go, go.getChilds(node->childs_id) + (a ^ 4), depth + 1, a, node_1, ray_source - vd3d{ 0.0, 0.0, node_1.z }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode1D(8);
				break;
			case 5:
				procSubtreeZ(tzm, tz1, go, go.getChilds(node->childs_id) + (a ^ 5), depth + 1, a, node_1, ray_source - vd3d{ node_1.x, 0.0, node_1.z }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode1D(8);
				break;
			case 6:
				procSubtreeZ(tzm, tz1, go, go.getChilds(node->childs_id) + (a ^ 6), depth + 1, a, node_1, ray_source - vd3d{ 0.0, node_1.y, node_1.z }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode1D(8);
				break;
			case 7:
				procSubtreeZ(tzm, tz1, go, go.getChilds(node->childs_id) + (a ^ 7), depth + 1, a, node_1, ray_source - vd3d{ node_1.x, node_1.y, node_1.z }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode1D(8);
				break;
			}
		} while (curr_node < 8);
	}

	void procSubtreeXY(double tx0, double ty0, double tx1, double ty1, GameObject& go, svo_model* node, uint8_t depth, uint8_t a, vd3d node_1, vd3d ray_source, ray_render_info& ray_info) {
		ray_info.max_depth = std::max(ray_info.max_depth, depth);
		ray_info.func_calls++;

		if (tx1 < 0.0 || ty1 < 0.0)
			return;
		double param = std::max(tx0, ty0);
		if (
			node->is_terminal() ||
			std::abs(node_1.sum() / param) < config.min_view_size
		) {
			if (go.materials[node->material_id].color.a != 255)
				return;
			ray_info.model = node;
			ray_info.depth = depth;
			// max
			if (tx0 > ty0) {
				ray_info.param = param;
				ray_info.axis_param = std::abs(ray_info.ray_direction.x);
				ray_info.axis = axis::X;
			}
			else {
				ray_info.param = param;
				ray_info.axis_param = std::abs(ray_info.ray_direction.y);
				ray_info.axis = axis::Y;
			}
			return;
		}

		double txm = 0.5 * (tx0 + tx1);
		double tym = 0.5 * (ty0 + ty1);

		node_1 *= 0.5;

		uint8_t curr_node = firstNodeXY(tx0, ty0, txm, tym, ray_source.z, node_1.z);

		do {
			switch (curr_node) {
			case 0:
				procSubtreeXY(tx0, ty0, txm, tym, go, go.getChilds(node->childs_id) + (a), depth + 1, a, node_1, ray_source, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode2D(txm, 1, tym, 2);
				break;
			case 1:
				procSubtreeXY(txm, ty0, tx1, tym, go, go.getChilds(node->childs_id) + (a ^ 1), depth + 1, a, node_1, ray_source - vd3d{ node_1.x, 0.0, 0.0 }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode2D(tx1, 8, tym, 3);
				break;
			case 2:
				procSubtreeXY(tx0, tym, txm, ty1, go, go.getChilds(node->childs_id) + (a ^ 2), depth + 1, a, node_1, ray_source - vd3d{ 0.0, node_1.y, 0.0 }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode2D(txm, 3, ty1, 8);
				break;
			case 3:
				procSubtreeXY(txm, tym, tx1, ty1, go, go.getChilds(node->childs_id) + (a ^ 3), depth + 1, a, node_1, ray_source - vd3d{ node_1.x, node_1.y, 0.0 }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode2D(tx1, 8, ty1, 8);
				break;
			case 4:
				procSubtreeXY(tx0, ty0, txm, tym, go, go.getChilds(node->childs_id) + (a ^ 4), depth + 1, a, node_1, ray_source - vd3d{ 0.0, 0.0, node_1.z }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode2D(txm, 5, tym, 6);
				break;
			case 5:
				procSubtreeXY(txm, ty0, tx1, tym, go, go.getChilds(node->childs_id) + (a ^ 5), depth + 1, a, node_1, ray_source - vd3d{ node_1.x, 0.0, node_1.z }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode2D(tx1, 8, tym, 7);
				break;
			case 6:
				procSubtreeXY(tx0, tym, txm, ty1, go, go.getChilds(node->childs_id) + (a ^ 6), depth + 1, a, node_1, ray_source - vd3d{ 0.0, node_1.y, node_1.z }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode2D(txm, 7, ty1, 8);
				break;
			case 7:
				procSubtreeXY(txm, tym, tx1, ty1, go, go.getChilds(node->childs_id) + (a ^ 7), depth + 1, a, node_1, ray_source - vd3d{ node_1.x, node_1.y, node_1.z }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode2D(tx1, 8, ty1, 8);
				break;
			}
		} while (curr_node < 8);
	}

	void procSubtreeXZ(double tx0, double tz0, double tx1, double tz1, GameObject& go, svo_model* node, uint8_t depth, uint8_t a, vd3d node_1, vd3d ray_source, ray_render_info& ray_info) {
		ray_info.max_depth = std::max(ray_info.max_depth, depth);
		ray_info.func_calls++;

		if (tx1 < 0.0 || tz1 < 0.0)
			return;
		double param = std::max(tx0, tz0);
		if (
			node->is_terminal() ||
			std::abs(node_1.sum() / param) < config.min_view_size
		) {
			if (go.materials[node->material_id].color.a != 255)
				return;
			ray_info.model = node;
			ray_info.depth = depth;
			// max
			if (tx0 > tz0) {
				ray_info.param = param;
				ray_info.axis_param = std::abs(ray_info.ray_direction.x);
				ray_info.axis = axis::X;
			}
			else {
				ray_info.param = param;
				ray_info.axis_param = std::abs(ray_info.ray_direction.z);
				ray_info.axis = axis::Z;
			}
			return;
		}

		double txm = 0.5 * (tx0 + tx1);
		double tzm = 0.5 * (tz0 + tz1);

		node_1 *= 0.5;

		uint8_t curr_node = firstNodeXZ(tx0, tz0, txm, tzm, ray_source.y, node_1.y);

		do {
			switch (curr_node) {
			case 0:
				procSubtreeXZ(tx0, tz0, txm, tzm, go, go.getChilds(node->childs_id) + (a), depth + 1, a, node_1, ray_source, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode2D(txm, 1, tzm, 4);
				break;
			case 1:
				procSubtreeXZ(txm, tz0, tx1, tzm, go, go.getChilds(node->childs_id) + (a ^ 1), depth + 1, a, node_1, ray_source - vd3d{ node_1.x, 0.0, 0.0 }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode2D(tx1, 8, tzm, 5);
				break;
			case 2:
				procSubtreeXZ(tx0, tz0, txm, tzm, go, go.getChilds(node->childs_id) + (a ^ 2), depth + 1, a, node_1, ray_source - vd3d{ 0.0, node_1.y, 0.0 }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode2D(txm, 3, tzm, 6);
				break;
			case 3:
				procSubtreeXZ(txm, tz0, tx1, tzm, go, go.getChilds(node->childs_id) + (a ^ 3), depth + 1, a, node_1, ray_source - vd3d{ node_1.x, node_1.y, 0.0 }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode2D(tx1, 8, tzm, 7);
				break;
			case 4:
				procSubtreeXZ(tx0, tzm, txm, tz1, go, go.getChilds(node->childs_id) + (a ^ 4), depth + 1, a, node_1, ray_source - vd3d{ 0.0, 0.0, node_1.z }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode2D(txm, 5, tz1, 8);
				break;
			case 5:
				procSubtreeXZ(txm, tzm, tx1, tz1, go, go.getChilds(node->childs_id) + (a ^ 5), depth + 1, a, node_1, ray_source - vd3d{ node_1.x, 0.0, node_1.z }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode2D(tx1, 8, tz1, 8);
				break;
			case 6:
				procSubtreeXZ(tx0, tzm, txm, tz1, go, go.getChilds(node->childs_id) + (a ^ 6), depth + 1, a, node_1, ray_source - vd3d{ 0.0, node_1.y, node_1.z }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode2D(txm, 7, tz1, 8);
				break;
			case 7:
				procSubtreeXZ(txm, tzm, tx1, tz1, go, go.getChilds(node->childs_id) + (a ^ 7), depth + 1, a, node_1, ray_source - vd3d{ node_1.x, node_1.y, node_1.z }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode2D(tx1, 8, tz1, 8);
				break;
			}
		} while (curr_node < 8);
	}

	void procSubtreeYZ(double ty0, double tz0, double ty1, double tz1, GameObject& go, svo_model* node, uint8_t depth, uint8_t a, vd3d node_1, vd3d ray_source, ray_render_info& ray_info) {
		ray_info.max_depth = std::max(ray_info.max_depth, depth);
		ray_info.func_calls++;

		if (ty1 < 0.0 || tz1 < 0.0)
			return;
		double param = std::max(ty0, tz0);
		if (
			node->is_terminal() ||
			std::abs(node_1.sum() / param) < config.min_view_size
		) {
			if (go.materials[node->material_id].color.a != 255)
				return;
			ray_info.model = node;
			ray_info.depth = depth;
			// max
			if (ty0 > tz0) {
				ray_info.param = param;
				ray_info.axis_param = std::abs(ray_info.ray_direction.y);
				ray_info.axis = axis::Y;
			}
			else {
				ray_info.param = param;
				ray_info.axis_param = std::abs(ray_info.ray_direction.z);
				ray_info.axis = axis::Z;
			}
			return;
		}

		double tym = 0.5 * (ty0 + ty1);
		double tzm = 0.5 * (tz0 + tz1);

		node_1 *= 0.5;

		uint8_t curr_node = firstNodeYZ(ty0, tz0, tym, tzm, ray_source.x, node_1.x);

		do {
			switch (curr_node) {
			case 0:
				procSubtreeYZ(ty0, tym, tz0, tzm, go, go.getChilds(node->childs_id) + (a), depth + 1, a, node_1, ray_source, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode2D(tym, 2, tzm, 4);
				break;
			case 1:
				procSubtreeYZ(ty0, tym, tz0, tzm, go, go.getChilds(node->childs_id) + (a ^ 1), depth + 1, a, node_1, ray_source - vd3d{ node_1.x, 0.0, 0.0 }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode2D(tym, 3, tzm, 5);
				break;
			case 2:
				procSubtreeYZ(tym, ty1, tz0, tzm, go, go.getChilds(node->childs_id) + (a ^ 2), depth + 1, a, node_1, ray_source - vd3d{ 0.0, node_1.y, 0.0 }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode2D(ty1, 8, tzm, 6);
				break;
			case 3:
				procSubtreeYZ(tym, ty1, tz0, tzm, go, go.getChilds(node->childs_id) + (a ^ 3), depth + 1, a, node_1, ray_source - vd3d{ node_1.x, node_1.y, 0.0 }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode2D(ty1, 8, tzm, 7);
				break;
			case 4:
				procSubtreeYZ(ty0, tym, tzm, tz1, go, go.getChilds(node->childs_id) + (a ^ 4), depth + 1, a, node_1, ray_source - vd3d{ 0.0, 0.0, node_1.z }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode2D(tym, 6, tz1, 8);
				break;
			case 5:
				procSubtreeYZ(ty0, tym, tzm, tz1, go, go.getChilds(node->childs_id) + (a ^ 5), depth + 1, a, node_1, ray_source - vd3d{ node_1.x, 0.0, node_1.z }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode2D(tym, 7, tz1, 8);
				break;
			case 6:
				procSubtreeYZ(tym, ty1, tzm, tz1, go, go.getChilds(node->childs_id) + (a ^ 6), depth + 1, a, node_1, ray_source - vd3d{ 0.0, node_1.y, node_1.z }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode2D(ty1, 8, tz1, 8);
				break;
			case 7:
				procSubtreeYZ(tym, ty1, tzm, tz1, go, go.getChilds(node->childs_id) + (a ^ 7), depth + 1, a, node_1, ray_source - vd3d{ node_1.x, node_1.y, node_1.z }, ray_info);
				if (ray_info.model != nullptr)
					return;
				curr_node = nextNode2D(ty1, 8, tz1, 8);
				break;
			}
		} while (curr_node < 8);
	}

	void procSubtreeXYZ(double tx0, double ty0, double tz0, double tx1, double ty1, double tz1, GameObject& go, svo_model* node, uint8_t depth, uint8_t a, vd3d node_1, ray_render_info& ray_info)
	{
		ray_info.max_depth = std::max(ray_info.max_depth, depth);
		ray_info.func_calls++;

		if (tx1 < 0.0 || ty1 < 0.0 || tz1 < 0.0)
			return;
		double param = std::max(tx0, std::max(ty0, tz0));
		if (
			node->is_terminal() ||
			std::abs(node_1.sum() / param) < config.min_view_size
		) {
			if (go.materials[node->material_id].color.a != 255)
				return;
			ray_info.model = node;
			ray_info.depth = depth;
			// max
			if (tx0 > ty0) {
				if (tx0 > tz0) {
					ray_info.param = param;
					ray_info.axis_param = std::abs(ray_info.ray_direction.x);
					ray_info.axis = axis::X;
				}
				else {
					ray_info.param = param;
					ray_info.axis_param = std::abs(ray_info.ray_direction.z);
					ray_info.axis = axis::Z;
				}
			}
			else if (ty0 > tz0) {
				ray_info.param = param;
				ray_info.axis_param = std::abs(ray_info.ray_direction.y);
				ray_info.axis = axis::Y;
			}
			else {
				ray_info.param = param;
				ray_info.axis_param = std::abs(ray_info.ray_direction.z);
				ray_info.axis = axis::Z;
			}
			return;
		}

		double txm = 0.5 * (tx0 + tx1);
		double tym = 0.5 * (ty0 + ty1);
		double tzm = 0.5 * (tz0 + tz1);

		uint8_t currNode = firstNodeXYZ(tx0, ty0, tz0, txm, tym, tzm);

		node_1 *= 0.5;

		do {
			switch (currNode)
			{
			case 0:
				procSubtreeXYZ(tx0, ty0, tz0, txm, tym, tzm, go, go.getChilds(node->childs_id) + (a), depth + 1, a, node_1, ray_info);
				if (ray_info.model != nullptr)
					return;
				currNode = nextNode3D(txm, 1, tym, 2, tzm, 4);
				break;
			case 1:
				procSubtreeXYZ(txm, ty0, tz0, tx1, tym, tzm, go, go.getChilds(node->childs_id) + (1 ^ a), depth + 1, a, node_1, ray_info);
				if (ray_info.model != nullptr)
					return;
				currNode = nextNode3D(tx1, 8, tym, 3, tzm, 5);
				break;
			case 2:
				procSubtreeXYZ(tx0, tym, tz0, txm, ty1, tzm, go, go.getChilds(node->childs_id) + (2 ^ a), depth + 1, a, node_1, ray_info);
				if (ray_info.model != nullptr)
					return;
				currNode = nextNode3D(txm, 3, ty1, 8, tzm, 6);
				break;
			case 3:
				procSubtreeXYZ(txm, tym, tz0, tx1, ty1, tzm, go, go.getChilds(node->childs_id) + (3 ^ a), depth + 1, a, node_1, ray_info);
				if (ray_info.model != nullptr)
					return;
				currNode = nextNode3D(tx1, 8, ty1, 8, tzm, 7);
				break;
			case 4:
				procSubtreeXYZ(tx0, ty0, tzm, txm, tym, tz1, go, go.getChilds(node->childs_id) + (4 ^ a), depth + 1, a, node_1, ray_info);
				if (ray_info.model != nullptr)
					return;
				currNode = nextNode3D(txm, 5, tym, 6, tz1, 8);
				break;
			case 5:
				procSubtreeXYZ(txm, ty0, tzm, tx1, tym, tz1, go, go.getChilds(node->childs_id) + (5 ^ a), depth + 1, a, node_1, ray_info);
				if (ray_info.model != nullptr)
					return;
				currNode = nextNode3D(tx1, 8, tym, 7, tz1, 8);
				break;
			case 6:
				procSubtreeXYZ(tx0, tym, tzm, txm, ty1, tz1, go, go.getChilds(node->childs_id) + (6 ^ a), depth + 1, a, node_1, ray_info);
				if (ray_info.model != nullptr)
					return;
				currNode = nextNode3D(txm, 7, ty1, 8, tz1, 8);
				break;
			case 7:
				procSubtreeXYZ(txm, tym, tzm, tx1, ty1, tz1, go, go.getChilds(node->childs_id) + (7 ^ a), depth + 1, a, node_1, ray_info);
				if (ray_info.model != nullptr)
					return;
				currNode = 8;
				break;
			}
		} while (currNode < 8);
	}

	inline uint8_t firstNodeX(vd3d& ray_source, vd3d& half_node) {
		uint8_t answer = 0;
		if (ray_source.y >= half_node.y) answer |= 2;
		if (ray_source.z >= half_node.z) answer |= 4;
		return answer;
	}

	inline uint8_t firstNodeY(vd3d& ray_source, vd3d& half_node) {
		uint8_t answer = 0;
		if (ray_source.x >= half_node.x) answer |= 1;
		if (ray_source.z >= half_node.z) answer |= 4;
		return answer;
	}

	inline uint8_t firstNodeZ(vd3d& ray_source, vd3d& half_node) {
		uint8_t answer = 0;
		if (ray_source.x >= half_node.x) answer |= 1;
		if (ray_source.y >= half_node.y) answer |= 2;
		return answer;
	}

	inline uint8_t firstNodeXY(double tx0, double ty0, double txm, double tym, double ray_source_z, double half_node_z) {
		uint8_t answer = 0;
		if (ray_source_z >= half_node_z) answer |= 4;
		if (tym < tx0) return (answer | 2);
		if (txm < ty0) return (answer | 1);
		return answer;
	}

	inline uint8_t firstNodeXZ(double tx0, double tz0, double txm, double tzm, double ray_source_y, double half_node_y) {
		uint8_t answer = 0;
		if (ray_source_y >= half_node_y) answer |= 2;
		if (tzm < tx0) return (answer | 4);
		if (txm < tz0) return (answer | 1);
		return answer;
	}

	inline uint8_t firstNodeYZ(double ty0, double tz0, double tym, double tzm, double ray_source_x, double half_node_x) {
		uint8_t answer = 0;
		if (ray_source_x >= half_node_x) answer |= 1;
		if (tzm < ty0) return (answer | 4);
		if (tym < tz0) return (answer | 2);
		return answer;
	}

	inline uint8_t firstNodeXYZ(double tx0, double ty0, double tz0, double txm, double tym, double tzm)
	{
		uint8_t answer = 0;
		if (tx0 > ty0)
		{
			if (tx0 > tz0)
			{
				if (tym < tx0)
					answer |= 2;
				if (tzm < tx0)
					answer |= 4;
				return answer;
			}
		}
		else
		{
			if (ty0 > tz0)
			{
				if (txm < ty0)
					answer |= 1;
				if (tzm < ty0)
					answer |= 4;
				return answer;
			}
		}
		if (txm < tz0)
			answer |= 1;
		if (tym < tz0)
			answer |= 2;
		return answer;
	}

	inline uint8_t nextNode1D(uint8_t x) {
		return x;
	}

	inline uint8_t nextNode2D(double tx, uint8_t x, double ty, uint8_t y) {
		return tx < ty ? x : y;
	}

	inline uint8_t nextNode3D(double txm, uint8_t x, double tym, uint8_t y, double tzm, uint8_t z)
	{
		if (txm < tym)
		{
			if (txm < tzm)
				return x;
		}
		else
		{
			if (tym < tzm)
				return y;
		}
		return z;
	}
	

protected:
	void setWorkerThreads(uint8_t count) {
		for (auto& w : workers) {
			w.alive = false;
		}
		work_cv.notify_all();
		for (auto& w : workers) {
			w.thread.join();
		}

		workers.clear();
		workers.resize(count);

		for (int i = 0; i < workers.size(); i++) {
			workers[i].thread = std::thread(&worker::job, &workers[i], std::ref(work_cv), std::ref(complete_cv));

			workers[i]._job = [&, height = ScreenHeight(), width = ScreenWidth(), start = i, count = workers.size()](vd3d luCorner, vd3d deltaX, vd3d deltaY) {
				for (uint32_t y = 0; y < height; y++)
				{
					vd3d left = luCorner + deltaY * y;
					for (uint32_t x = start; x < width; x += count)
					{
						vd3d rayPoint = left + deltaX * x;
						vd3d rayDir = rayPoint - player_view.pos;

#ifdef SVO_RENDER_V2
						ray_render_info ret = rayParameterV2(testObject, rayPoint, rayDir);
#else
						ray_render_info ret = rayParameter(testObject, rayPoint, rayDir);
#endif
						std::invoke(*render_func_it, *this, x, y, std::ref(testObject), std::ref(ret));

						if (x == stop_x && y == stop_y) {
							stop_x = ScreenWidth();
							stop_y = ScreenHeight();
							//std::cout << "Pause" << std::endl;
							std::cout << ret << std::endl;
						}
					}
				}

				work_in_progress--;
			};
		}

		for (auto& w : workers) {
			while (!w.started) { std::this_thread::yield(); }
			{ std::unique_lock<std::mutex> ul(w.mut); }
		}
	}


protected:
	void RenderDefault(uint32_t x, uint32_t y, GameObject& go, ray_render_info& rri) {
		if (rri.model)
			Draw(x, y, go.materials[rri.model->material_id].color);
		else
			Draw(x, y, olc::BLANK);
	}

	void RenderLength(uint32_t x, uint32_t y, GameObject& go, ray_render_info& rri) {
		float n = rri.param;
		float a = 0.0001f;
		if (rri.model)
			Draw(x, y, olc::PixelF(0.5f * sinf(a * n) + 0.5f, 0.5f * sin(a * n + 2.094f) + 0.5f, 0.5f * sinf(a * n + 4.188f) + 0.5f));
		else
			Draw(x, y, olc::BLANK);
	}

	void RenderDepth(uint32_t x, uint32_t y, GameObject& go, ray_render_info& rri) {
		float n = rri.depth;
		float a = 0.1f;
		if (rri.model)
			Draw(x, y, olc::PixelF(0.5f * sinf(a * n) + 0.5f, 0.5f * sin(a * n + 2.094f) + 0.5f, 0.5f * sinf(a * n + 4.188f) + 0.5f));
		else
			Draw(x, y, olc::BLANK);
	}

	void RenderMaxDepth(uint32_t x, uint32_t y, GameObject& go, ray_render_info& rri) {
		float n = rri.max_depth;
		float a = 0.1f;
		Draw(x, y, olc::PixelF(0.5f * sinf(a * n) + 0.5f, 0.5f * sin(a * n + 2.094f) + 0.5f, 0.5f * sinf(a * n + 4.188f) + 0.5f));
	}

	void RenderFuncCalls(uint32_t x, uint32_t y, GameObject& go, ray_render_info& rri) {
		float n = rri.func_calls;
		float a = 0.1f;
		Draw(x, y, olc::PixelF(0.5f * sinf(a * n) + 0.5f, 0.5f * sin(a * n + 2.094f) + 0.5f, 0.5f * sinf(a * n + 4.188f) + 0.5f));
	}
};

int main(int argc, char *argv[])
{
	SVOGame game;
	if (game.Construct(600, 400, 2, 2))
		game.Start();
	return 0;
}
