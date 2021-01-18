#define OLC_PGE_APPLICATION
#define _USE_MATH_DEFINES
#include "olcPixelGameEngine.h"


//struct svo_base {
//
//};


//struct svo_object {
//	uint8_t size;
//	olc::Pixel color;
//	std::vector<svo_object> childs;
//	svo_object* parent;
//	uint8_t used_childs_mask;
//
//	svo_object* create_childs() {
//		childs.resize(8);
//		childs[0].parent = this;
//		childs[1].parent = this;
//		childs[2].parent = this;
//		childs[3].parent = this;
//		childs[4].parent = this;
//		childs[5].parent = this;
//		childs[6].parent = this;
//		childs[7].parent = this;
//	}
//};
//
//
//class game_object {
//	uint32_t xSize, ySize, zSize;
//	double x, y, z;
//};
//
//
//template<typename T>
//struct octree_node {
//	std::vector<T> items;
//	std::vector<octree_node> childs;
//	octree_node* parent = nullptr;
//
//	octree_node* create_childs() {
//		childs.resize(8);
//		childs[0].parent = this;
//		childs[1].parent = this;
//		childs[2].parent = this;
//		childs[3].parent = this;
//		childs[4].parent = this;
//		childs[5].parent = this;
//		childs[6].parent = this;
//		childs[7].parent = this;
//	}
//};
//
//
//class Space3D {
//
//};


class camera {
	double x, y, z;
	double hAngle, vAngle;
	const double move_speed = 0.1;
	const double turn_speed = 0.01;

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

	double X() {
		return x;
	}
	double Y() {
		return y;
	}
	double Z() {
		return z;
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

		std::stringstream coords;
		coords << "X: " << player_view.X() << '\n';
		coords << "Y: " << player_view.Y() << '\n';
		coords << "Z: " << player_view.Z() << '\n';
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