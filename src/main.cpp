#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "GLWrapper.h"
#include "SceneManager.h"
#include "Surface.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <ctime>

static int wind_width = 1280;
static int wind_height = 720;

void update_scene(scene_container& scene, float delta, float time);

// 按照三角形中心排序 -- 比较函数
bool cmpx(const rt_triangle& t1, const rt_triangle& t2) {
	glm::vec3 center1 = (t1.p1 + t1.p2 + t1.p3) / glm::vec3(3, 3, 3);
	glm::vec3 center2 = (t2.p1 + t2.p2 + t2.p3) / glm::vec3(3, 3, 3);
	return center1.x < center2.x;
}
bool cmpy(const rt_triangle& t1, const rt_triangle& t2) {
	glm::vec3 center1 = (t1.p1 + t1.p2 + t1.p3) / glm::vec3(3, 3, 3);
	glm::vec3 center2 = (t2.p1 + t2.p2 + t2.p3) / glm::vec3(3, 3, 3);
	return center1.y < center2.y;
}
bool cmpz(const rt_triangle& t1, const rt_triangle& t2) {
	glm::vec3 center1 = (t1.p1 + t1.p2 + t1.p3) / glm::vec3(3, 3, 3);
	glm::vec3 center2 = (t2.p1 + t2.p2 + t2.p3) / glm::vec3(3, 3, 3);
	return center1.z < center2.z;
}

void readObj(std::string filepath, std::vector<rt_triangle>& triangles, rt_material material, glm::mat4 trans, bool smoothNormal) {

	// 顶点位置，索引
	std::vector<glm::vec3> vertices;
	std::vector<GLuint> indices;

	// 打开文件流
	std::ifstream fin(filepath);
	std::string line;
	if (!fin.is_open()) {
		std::cout << "文件 " << filepath << " 打开失败" << std::endl;
		exit(-1);
	}

	// 计算 AABB 盒，归一化模型大小
	float maxx = -11451419.19;
	float maxy = -11451419.19;
	float maxz = -11451419.19;
	float minx = 11451419.19;
	float miny = 11451419.19;
	float minz = 11451419.19;

	// 按行读取
	while (std::getline(fin, line)) {
		std::istringstream sin(line);   // 以一行的数据作为 string stream 解析并且读取
		std::string type;
		GLfloat x, y, z;
		int v0, v1, v2;
		int vn0, vn1, vn2;
		int vt0, vt1, vt2;
		char slash;

		// 统计斜杆数目，用不同格式读取
		int slashCnt = 0;
		for (int i = 0; i < line.length(); i++) {
			if (line[i] == '/') slashCnt++;
		}

		// 读取obj文件
		sin >> type;
		if (type == "v") {
			sin >> x >> y >> z;
			vertices.push_back(glm::vec3(x, y, z));
			maxx = std::max(maxx, x); maxy = std::max(maxx, y); maxz = std::max(maxx, z);
			minx = std::min(minx, x); miny = std::min(minx, y); minz = std::min(minx, z);
		}
		if (type == "f") {
			if (slashCnt == 6) {
				sin >> v0 >> slash >> vt0 >> slash >> vn0;
				sin >> v1 >> slash >> vt1 >> slash >> vn1;
				sin >> v2 >> slash >> vt2 >> slash >> vn2;
			}
			else if (slashCnt == 3) {
				sin >> v0 >> slash >> vt0;
				sin >> v1 >> slash >> vt1;
				sin >> v2 >> slash >> vt2;
			}
			else {
				sin >> v0 >> v1 >> v2;
			}
			indices.push_back(v0 - 1);
			indices.push_back(v1 - 1);
			indices.push_back(v2 - 1);
		}
	}

	// 模型大小归一化
	float lenx = maxx - minx;
	float leny = maxy - miny;
	float lenz = maxz - minz;
	float maxaxis = std::max(lenx, std::max(leny, lenz));
	for (auto& v : vertices) {
		v.x /= maxaxis;
		v.y /= maxaxis;
		v.z /= maxaxis;
	}

	// 通过矩阵进行坐标变换
	for (auto& v : vertices) {
		glm::vec4 vv = glm::vec4(v.x, v.y, v.z, 1);
		vv = trans * vv;
		v = glm::vec3(vv.x, vv.y, vv.z);
	}

	// 生成法线
	std::vector<glm::vec3> normals(vertices.size(), glm::vec3(0, 0, 0));
	for (int i = 0; i < indices.size(); i += 3) {
		glm::vec3 p1 = vertices[indices[i]];
		glm::vec3 p2 = vertices[indices[i + 1]];
		glm::vec3 p3 = vertices[indices[i + 2]];
		glm::vec3 n = glm::normalize(glm::cross(p2 - p1, p3 - p1));
		normals[indices[i]] += n;
		normals[indices[i + 1]] += n;
		normals[indices[i + 2]] += n;
	}

	// 构建 rt_triangle 对象数组
	int offset = triangles.size();  // 增量更新
	triangles.resize(offset + indices.size() / 3);
	for (int i = 0; i < indices.size(); i += 3) {
		rt_triangle& t = triangles[offset + i / 3];
		// 传顶点属性
		t.p1 = vertices[indices[i]];
		t.p2 = vertices[indices[i + 1]];
		t.p3 = vertices[indices[i + 2]];
		if (!smoothNormal) {
			glm::vec3 n = glm::normalize(glm::cross(t.p2 - t.p1, t.p3 - t.p1));
			t.n1 = n; t.n2 = n; t.n3 = n;
		}
		else {
			t.n1 = normalize(normals[indices[i]]);
			t.n2 = normalize(normals[indices[i + 1]]);
			t.n3 = normalize(normals[indices[i + 2]]);
		}

		// 传材质
		t.material = material;
	}
}

// SAH 优化构建 BVH
int buildBVHwithSAH(std::vector<rt_triangle>& triangles, std::vector<rt_BVHNode>& nodes, int l, int r, int n) {
	using namespace glm;
	if (l > r) return 0;

	nodes.push_back(rt_BVHNode());
	int id = nodes.size() - 1;
	nodes[id].left = nodes[id].right = nodes[id].n = nodes[id].index = 0;
	nodes[id].AA = vec3(1145141919, 1145141919, 1145141919);
	nodes[id].BB = vec3(-1145141919, -1145141919, -1145141919);

	// 计算 AABB
	for (int i = l; i <= r; i++) {
		// 最小点 AA
		float minx = min(triangles[i].p1.x, min(triangles[i].p2.x, triangles[i].p3.x));
		float miny = min(triangles[i].p1.y, min(triangles[i].p2.y, triangles[i].p3.y));
		float minz = min(triangles[i].p1.z, min(triangles[i].p2.z, triangles[i].p3.z));
		nodes[id].AA.x = min(nodes[id].AA.x, minx);
		nodes[id].AA.y = min(nodes[id].AA.y, miny);
		nodes[id].AA.z = min(nodes[id].AA.z, minz);
		// 最大点 BB
		float maxx = max(triangles[i].p1.x, max(triangles[i].p2.x, triangles[i].p3.x));
		float maxy = max(triangles[i].p1.y, max(triangles[i].p2.y, triangles[i].p3.y));
		float maxz = max(triangles[i].p1.z, max(triangles[i].p2.z, triangles[i].p3.z));
		nodes[id].BB.x = max(nodes[id].BB.x, maxx);
		nodes[id].BB.y = max(nodes[id].BB.y, maxy);
		nodes[id].BB.z = max(nodes[id].BB.z, maxz);
	}

	// 不多于 n 个三角形 返回叶子节点
	if ((r - l + 1) <= n) {
		nodes[id].n = r - l + 1;
		nodes[id].index = l;
		return id;
	}

	// 否则递归建树
	#define INF 10000000.0
	float Cost = INF;
	int Axis = 0;
	int Split = (l + r) / 2;
	for (int axis = 0; axis < 3; axis++) {
		// 分别按 x，y，z 轴排序
		if (axis == 0) std::sort(&triangles[0] + l, &triangles[0] + r + 1, cmpx);
		if (axis == 1) std::sort(&triangles[0] + l, &triangles[0] + r + 1, cmpy);
		if (axis == 2) std::sort(&triangles[0] + l, &triangles[0] + r + 1, cmpz);

		// leftMax[i]: [l, i] 中最大的 xyz 值
		// leftMin[i]: [l, i] 中最小的 xyz 值
		std::vector<vec3> leftMax(r - l + 1, vec3(-INF, -INF, -INF));
		std::vector<vec3> leftMin(r - l + 1, vec3(INF, INF, INF));
		// 计算前缀 注意 i-l 以对齐到下标 0
		for (int i = l; i <= r; i++) {
			rt_triangle& t = triangles[i];
			int bias = (i == l) ? 0 : 1;  // 第一个元素特殊处理

			leftMax[i - l].x = max(leftMax[i - l - bias].x, max(t.p1.x, max(t.p2.x, t.p3.x)));
			leftMax[i - l].y = max(leftMax[i - l - bias].y, max(t.p1.y, max(t.p2.y, t.p3.y)));
			leftMax[i - l].z = max(leftMax[i - l - bias].z, max(t.p1.z, max(t.p2.z, t.p3.z)));

			leftMin[i - l].x = min(leftMin[i - l - bias].x, min(t.p1.x, min(t.p2.x, t.p3.x)));
			leftMin[i - l].y = min(leftMin[i - l - bias].y, min(t.p1.y, min(t.p2.y, t.p3.y)));
			leftMin[i - l].z = min(leftMin[i - l - bias].z, min(t.p1.z, min(t.p2.z, t.p3.z)));
		}

		// rightMax[i]: [i, r] 中最大的 xyz 值
		// rightMin[i]: [i, r] 中最小的 xyz 值
		std::vector<vec3> rightMax(r - l + 1, vec3(-INF, -INF, -INF));
		std::vector<vec3> rightMin(r - l + 1, vec3(INF, INF, INF));
		// 计算后缀 注意 i-l 以对齐到下标 0
		for (int i = r; i >= l; i--) {
			rt_triangle& t = triangles[i];
			int bias = (i == r) ? 0 : 1;  // 第一个元素特殊处理

			rightMax[i - l].x = max(rightMax[i - l + bias].x, max(t.p1.x, max(t.p2.x, t.p3.x)));
			rightMax[i - l].y = max(rightMax[i - l + bias].y, max(t.p1.y, max(t.p2.y, t.p3.y)));
			rightMax[i - l].z = max(rightMax[i - l + bias].z, max(t.p1.z, max(t.p2.z, t.p3.z)));

			rightMin[i - l].x = min(rightMin[i - l + bias].x, min(t.p1.x, min(t.p2.x, t.p3.x)));
			rightMin[i - l].y = min(rightMin[i - l + bias].y, min(t.p1.y, min(t.p2.y, t.p3.y)));
			rightMin[i - l].z = min(rightMin[i - l + bias].z, min(t.p1.z, min(t.p2.z, t.p3.z)));
		}

		// 遍历寻找分割
		float cost = INF;
		int split = l;
		for (int i = l; i <= r - 1; i++) {
			float lenx, leny, lenz;
			// 左侧 [l, i]
			vec3 leftAA = leftMin[i - l];
			vec3 leftBB = leftMax[i - l];
			lenx = leftBB.x - leftAA.x;
			leny = leftBB.y - leftAA.y;
			lenz = leftBB.z - leftAA.z;
			float leftS = 2.0 * ((lenx * leny) + (lenx * lenz) + (leny * lenz));
			float leftCost = leftS * (i - l + 1);

			// 右侧 [i+1, r]
			vec3 rightAA = rightMin[i + 1 - l];
			vec3 rightBB = rightMax[i + 1 - l];
			lenx = rightBB.x - rightAA.x;
			leny = rightBB.y - rightAA.y;
			lenz = rightBB.z - rightAA.z;
			float rightS = 2.0 * ((lenx * leny) + (lenx * lenz) + (leny * lenz));
			float rightCost = rightS * (r - i);

			// 记录每个分割的最小答案
			float totalCost = leftCost + rightCost;
			if (totalCost < cost) {
				cost = totalCost;
				split = i;
			}
		}
		// 记录每个轴的最佳答案
		if (cost < Cost) {
			Cost = cost;
			Axis = axis;
			Split = split;
		}
	}

	// 按最佳轴分割
	if (Axis == 0) std::sort(&triangles[0] + l, &triangles[0] + r + 1, cmpx);
	if (Axis == 1) std::sort(&triangles[0] + l, &triangles[0] + r + 1, cmpy);
	if (Axis == 2) std::sort(&triangles[0] + l, &triangles[0] + r + 1, cmpz);

	// 递归
	int left = buildBVHwithSAH(triangles, nodes, l, Split, n);
	int right = buildBVHwithSAH(triangles, nodes, Split + 1, r, n);

	nodes[id].left = left;
	nodes[id].right = right;

	return id;
}

glm::mat4 getTransformMatrix(glm::vec3 rotateCtrl, glm::vec3 translateCtrl, glm::vec3 scaleCtrl) {
	using namespace glm;
	glm::mat4 unit(    // 单位矩阵
		glm::vec4(1, 0, 0, 0),
		glm::vec4(0, 1, 0, 0),
		glm::vec4(0, 0, 1, 0),
		glm::vec4(0, 0, 0, 1)
	);
	mat4 scale = glm::scale(unit, scaleCtrl);
	mat4 translate = glm::translate(unit, translateCtrl);
	mat4 rotate = unit;
	rotate = glm::rotate(rotate, glm::radians(rotateCtrl.x), glm::vec3(1, 0, 0));
	rotate = glm::rotate(rotate, glm::radians(rotateCtrl.y), glm::vec3(0, 1, 0));
	rotate = glm::rotate(rotate, glm::radians(rotateCtrl.z), glm::vec3(0, 0, 1));

	mat4 model = translate * rotate * scale;
	return model;
}

namespace update {
	int jupiter = -1,
		saturn = -1,
		saturn_rings = -1,
		mars = -1,
		box = -1,
		torus = -1;
}

const glm::quat saturn_pitch = glm::quat(glm::vec3(glm::radians(15.f), 0, 0));

int main()
{
	GLWrapper glWrapper(wind_width, wind_height, false);
	// fullscreen
	//GLWrapper glWrapper(true);
	// window with monitor resolution
	//GLWrapper glWrapper(false);

	// set SMAA quality preset
	glWrapper.enable_SMAA(ULTRA);
	
	glWrapper.init_window();
	glfwSwapInterval(1); // vsync
	wind_width = glWrapper.getWidth();
	wind_height = glWrapper.getHeight();

	// fix ray direction issues
	if (wind_width % 2 == 1) wind_width++;
	if (wind_height % 2 == 1) wind_height++;

	scene_container scene = {};

	scene.scene = SceneManager::create_scene(wind_width, wind_height);
	scene.scene.camera_pos = { 0, 0, -5 };
	scene.shadow_ambient = glm::vec3{ 0.1, 0.1, 0.1 };
	scene.ambient_color = glm::vec3{ 0.025, 0.025, 0.025 };

	// lights
	scene.lights_point.push_back(SceneManager::create_light_point({ 3, 5, 0, 0.1 }, { 1, 1, 1 }, 25.5));
	scene.lights_direct.push_back(SceneManager::create_light_direct({ 3, -1, 1 }, { 1, 1, 1 }, 1.5));

	// blue sphere
	scene.spheres.push_back(SceneManager::create_sphere({ 2, 0, 6 }, 1,
		SceneManager::create_material({ 0, 0, 1 }, 50, 0.35)));
	// red sphere
	scene.spheres.push_back(SceneManager::create_sphere({ -1, 0, 6 }, 1,
		SceneManager::create_material({ 1, 0, 0 }, 100, 0.1), true));
	// transparent sphere
	scene.spheres.push_back(SceneManager::create_sphere({ 0.5, 2, 6 }, 1,
		SceneManager::create_material({ 1, 1, 1 }, 200, 0.1, 1.125, { 0, 0, 0 }, 1), true));

	std::vector<rt_triangle> triangles;
	readObj("cube.obj", triangles, SceneManager::create_material({ 1, 0, 0 }, 100, 0.1), getTransformMatrix(glm::vec3(0, 0, 0), glm::vec3(0.3, -1.6, 0), glm::vec3(1.5, 1.5, 1.5)), true);
	scene.triangles = triangles;
	int nTriangles = triangles.size();
	std::cout << "模型读取完成: 共 " << nTriangles << " 个三角形" << std::endl;
	
	rt_BVHNode testnode;
	testnode.left = 255;
	testnode.right = 128;
	testnode.n = 30;
	testnode.AA = glm::vec3(1, 1, 0);
	testnode.BB = glm::vec3(0, 1, 0);
	std::vector<rt_BVHNode> nodes{ testnode };
	//buildbvh(triangles, nodes, 0, triangles.size() - 1, 8);
	buildBVHwithSAH(triangles, nodes, 0, triangles.size() - 1, 8);
	scene.BVHNodes = nodes;
	int nnodes = nodes.size();
	std::cout << "bvh 建立完成: 共 " << nnodes << " 个节点" << std::endl;

	 //jupiter
	 //jupiter 的位置定义在update_scene中
	rt_sphere jupiter = SceneManager::create_sphere({}, 5000,
		SceneManager::create_material({}, 0, 0.0f));
	jupiter.textureNum = 1;
	scene.spheres.push_back(jupiter);
	update::jupiter = scene.spheres.size() - 1;

	// saturn
	const int saturnRadius = 4150;
	rt_sphere saturn = SceneManager::create_sphere({}, saturnRadius,
		SceneManager::create_material({}, 0, 0.0f));
	saturn.textureNum = 2;
	saturn.quat_rotation = saturn_pitch;
	scene.spheres.push_back(saturn);
	update::saturn = scene.spheres.size() - 1;

	// mars
	rt_sphere mars = SceneManager::create_sphere({}, 500,
		SceneManager::create_material({}, 0, 0.0f));
	mars.textureNum = 3;
	scene.spheres.push_back(mars);
	update::mars = scene.spheres.size() - 1;

	// ring
	{
		rt_ring ring = SceneManager::create_ring({}, saturnRadius * 1.1166, saturnRadius * 2.35,
			SceneManager::create_material({}, 0, 0));
		ring.textureNum = 4;
		ring.quat_rotation = glm::angleAxis(glm::radians(90.f), glm::vec3(1, 0, 0)) * saturn_pitch;
		scene.rings.push_back(ring);
		update::saturn_rings = scene.rings.size() - 1;
	}

	// floor
	scene.boxes.push_back(SceneManager::create_box({ 0, -1.2, 6 }, { 10, 0.2, 5 },
		SceneManager::create_material({ 1, 0.6, 0 }, 100, 0.05)));
	// box
	rt_box box = SceneManager::create_box({ 8, 1, 6 }, { 1, 1, 1 },
		SceneManager::create_material({ 0.8,0.7,0 }, 50, 0.0));
	box.textureNum = 5;
	scene.boxes.push_back(box);
	update::box = scene.boxes.size() - 1;

	// *** beware! torus calculations is the heaviest part of rendering
	// *** comment next line if you have performance issues
	// torus
	rt_torus torus = SceneManager::create_torus({ -9, 0.5, 6 }, { 1.0, 0.5 },
		SceneManager::create_material({ 0.5, 0.4, 1 }, 200, 0.2));
	torus.quat_rotation = glm::quat(glm::vec3(glm::radians(45.f), 0, 0));
	scene.toruses.push_back(torus);
	update::torus = scene.toruses.size() - 1;

	// cone
	rt_material coneMaterial = SceneManager::create_material({ 234 / 255.0f, 17 / 255.0f, 82 / 255.0f }, 200, 0.2);
	rt_surface cone = SurfaceFactory::GetEllipticCone(1 / 3.0f, 1 / 3.0f, 1, coneMaterial);
	cone.pos = { -5, 4, 6 };
	cone.quat_rotation = glm::quat(glm::vec3(glm::radians(90.f), 0, 0));
	cone.yMin = -1;
	cone.yMax = 4;
	scene.surfaces.push_back(cone);

	// cylinder
	rt_material cylinderMaterial = SceneManager::create_material({ 200 / 255.0f, 255 / 255.0f, 0 / 255.0f }, 200, 0.2);
	rt_surface cylinder = SurfaceFactory::GetEllipticCylinder(1 / 2.0f, 1 / 2.0f, cylinderMaterial);
	cylinder.pos = { 5, 0, 6 };
	cylinder.quat_rotation = glm::quat(glm::vec3(glm::radians(90.f), 0, 0));
	cylinder.yMin = -1;
	cylinder.yMax = 1;
	scene.surfaces.push_back(cylinder);

	rt_defines defines = scene.get_defines();
	glWrapper.init_shaders(defines);

	//std::vector<std::string> faces =
	//{
	//	ASSETS_DIR "/textures/sb_nebula/GalaxyTex_PositiveX.jpg",
	//	ASSETS_DIR "/textures/sb_nebula/GalaxyTex_NegativeX.jpg",
	//	ASSETS_DIR "/textures/sb_nebula/GalaxyTex_PositiveY.jpg",
	//	ASSETS_DIR "/textures/sb_nebula/GalaxyTex_NegativeY.jpg",
	//	ASSETS_DIR "/textures/sb_nebula/GalaxyTex_PositiveZ.jpg",
	//	ASSETS_DIR "/textures/sb_nebula/GalaxyTex_NegativeZ.jpg"
	//};
	std::vector<std::string> faces =
	{
		ASSETS_DIR "/textures/skybox1/px.jpg",
		ASSETS_DIR "/textures/skybox1/nx.jpg",
		ASSETS_DIR "/textures/skybox1/py.jpg",
		ASSETS_DIR "/textures/skybox1/ny.jpg",
		ASSETS_DIR "/textures/skybox1/pz.jpg",
		ASSETS_DIR "/textures/skybox1/nz.jpg"
	};

	glWrapper.set_skybox(GLWrapper::load_cubemap(faces, false));
	//glWrapper.set_skybox(GLWrapper::load_cubemap(faces, true));

	auto jupiterTex = glWrapper.load_texture(1, "8k_jupiter.jpg", "texture_sphere_1");
	auto saturnTex = glWrapper.load_texture(2, "8k_saturn.jpg", "texture_sphere_2");
	auto marsTex = glWrapper.load_texture(3, "2k_mars.jpg", "texture_sphere_3");
	auto ringTex = glWrapper.load_texture(4, "8k_saturn_ring_alpha.png", "texture_ring");
	auto boxTex = glWrapper.load_texture(5, "xiaohui.png", "texture_box");

	SceneManager scene_manager(wind_width, wind_height, &scene, &glWrapper);
	scene_manager.init();

	float currentTime = static_cast<float>(glfwGetTime());
	float lastFramesPrint = currentTime;
	float framesCount = 0;

	while (!glfwWindowShouldClose(glWrapper.window))
	{
		framesCount++;
		float newTime = static_cast<float>(glfwGetTime());
		float deltaTime = newTime - currentTime;
		currentTime = newTime;

		if (newTime - lastFramesPrint > 1.0f)
		{
			std::cout << "FPS: " << framesCount << std::endl;
			lastFramesPrint = newTime;
			framesCount = 0;
		}

		update_scene(scene, deltaTime, newTime);
		scene_manager.update(deltaTime);
		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, jupiterTex);
		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, saturnTex);
		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_2D, marsTex);
		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_2D, ringTex);
		glActiveTexture(GL_TEXTURE5);
		glBindTexture(GL_TEXTURE_2D, boxTex);
		glWrapper.draw();
		glfwSwapBuffers(glWrapper.window);
		glfwPollEvents();
	}

	glWrapper.stop(); // stop glfw, close window
	return 0;
}

void update_scene(scene_container& scene, float deltaTime, float time)
{
	if (update::jupiter != -1) {
		rt_sphere* jupiter = &scene.spheres[update::jupiter];
		const float jupiterSpeed = 0.02;
		jupiter->obj.x = cos(time * jupiterSpeed) * 20000;
		jupiter->obj.z = sin(time * jupiterSpeed) * 20000;

		jupiter->quat_rotation *= glm::angleAxis(deltaTime / 15, glm::vec3(0, 1, 0));
	}

	if (update::saturn != -1 && update::saturn_rings != -1) {
		rt_sphere* saturn = &scene.spheres[update::saturn];
		rt_ring* ring = &scene.rings[update::saturn_rings];
		const float speed = 0.0082;
		const float dist = 35000;
		const float offset = 1;

		saturn->obj.x = cos(time * speed + offset) * dist;
		saturn->obj.z = sin(time * speed + offset) * dist;

		glm::vec3 axis = glm::vec3(0, 1, 0) * saturn_pitch;
		saturn->quat_rotation *= glm::angleAxis(deltaTime / 10, axis);

		ring->pos.x = cos(time * speed + offset) * dist;
		ring->pos.z = sin(time * speed + offset) * dist;
	}

	if (update::mars != -1) {
		rt_sphere* mars = &scene.spheres[update::mars];
		const float marsSpeed = 0.05;
		mars->obj.x = cos(time * marsSpeed + 0.5f) * 10000;
		mars->obj.z = sin(time * marsSpeed + 0.5f) * 10000;
		mars->obj.y = -cos(time * marsSpeed) * 3000;
		mars->quat_rotation *= glm::angleAxis(deltaTime / 5, glm::vec3(0, 1, 0));
	}

	if (update::box != -1)
	{
		rt_box* box = &scene.boxes[update::box];
		glm::quat q = glm::angleAxis(deltaTime, glm::vec3(0.5774, 0.5774, 0.5774));
		box->quat_rotation *= q;
	}

	if (update::torus != -1)
	{
		rt_torus* torus = &scene.toruses[update::torus];
		torus->quat_rotation *= glm::angleAxis(deltaTime, glm::vec3(0, 1, 0));
	}
}
