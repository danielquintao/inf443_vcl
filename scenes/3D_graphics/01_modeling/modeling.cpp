
#include "modeling.hpp"
#include <iostream>

#ifdef SCENE_3D_GRAPHICS

// Add vcl namespace within the current one - Allows to use function from vcl library without explicitely preceeding their name with vcl::
using namespace vcl;




float evaluate_terrain_z(float u, float v);
vec3 evaluate_terrain(float u, float v);
mesh create_terrain();
mesh create_pyramid(float radius, float height, float z_offset, float rot);
mesh create_tronc(float radius, float height);
mesh create_foliage(float radius, float height, float z_offset);
hierarchy_mesh_drawable create_tree();
mesh mesh_skybox();
float neck_position(float t, float& t_max);
float leg_position(float t, float& t_max);
float camel_position(float t, float& t_max);



/** This function is called before the beginning of the animation loop
    It is used to initialize all part-specific data */
void scene_model::setup_data(std::map<std::string,GLuint>& shaders , scene_structure& scene, gui_structure& )
{
    // Create visual terrain surface
    terrain = create_terrain();
    terrain.uniform.color = { 241 / 255.0f,175 / 255.0f,0.0f };
    terrain.uniform.shading.specular = 0.0f; // non-specular terrain material

    // sand texture
    // Load a texture image on GPU and stores its ID
    terrain.texture_id = create_texture_gpu(image_load_png("scenes/3D_graphics/01_modeling/assets/sand.png"));
    
    // Setup initial camera mode and position
    scene.camera.camera_type = camera_control_spherical_coordinates;
    scene.camera.scale = 10.0f;
    scene.camera.apply_rotation(0.8,0,0,1.5f);

    // oasis small lake :
    water = mesh_drawable(mesh_primitive_quad({ -3.f, -3.f, 0 }, { -3.f, 3.f, 0 }, { 3.f, 3.f, 0 }, { 3.f, -3.f, 0 }));
    water.uniform.color = { 0.0f, 0.5f, 1.0f };
    water.uniform.color_alpha = 0.5f; // transparency
    water.uniform.shading.specular = 0.9f;

    // Pyramid :
    // GEOMETRIC PYRAMID:
    //pyramid = mesh_drawable(create_pyramid(6.0f, 4.0f, 0.0f, 0.25f));
    //pyramid.uniform.transform.translation = { 5,5,0 };
    // IMPORTED PYRAMID:
    pyramid = mesh_drawable(mesh_load_file_obj("scenes/3D_graphics/01_modeling/assets/super_pyramid.obj"));
    vcl::mat3 R_pyramid = vcl::rotation_from_axis_angle_mat3({ 1,0,0 }, 3.14f / 2); // only in case of imported object
    pyramid.uniform.transform.rotation = R_pyramid; // only in case of imported object
    pyramid.uniform.transform.scaling = 1.5; // only in case of imported object
    pyramid.uniform.transform.translation = { 5,5,1.4f }; // only in case of imported object
    // GENERAL PYRAMID SETTINGS:
    pyramid.uniform.color = {0.86f, 0.65f, 0.30f};
    pyramid.uniform.shading.specular = 0.0f;
    pyramid.uniform.shading.ambiant = 0.5f;
    

    // Camel
    // main part of camel:
    mesh_drawable camel_trunk = mesh_drawable(mesh_load_file_obj("scenes/3D_graphics/01_modeling/assets/camelo_trunk.obj"));
    mesh_drawable camel_head = mesh_drawable(mesh_load_file_obj("scenes/3D_graphics/01_modeling/assets/camelo_head.obj"));
    camel_trunk.uniform.color = { 1.0f, 0.75f, 0.2f };
    camel_trunk.uniform.shading.specular = 0.01f;
    camel_head.uniform.color = { 1.0f, 0.75f, 0.2f };
    camel_head.uniform.shading.specular = 0.01f;
    vcl::mat3 R_camel = vcl::rotation_from_axis_angle_mat3({ 1,0,0 }, 3.14f / 2);
    camel_trunk.uniform.transform.rotation = R_camel;
    camel_head.uniform.transform.rotation = R_camel;
    camel_head.uniform.transform.translation = { 0, 0.17f, 0 };
    // leg:
    mesh_drawable articulation_point = mesh_primitive_sphere(0.02f, { 0,0,0 }, 10, 10);
    mesh_drawable thigh = mesh_primitive_cylinder(0.025f, { 0,0,-0.2f }, { 0,0,0 }, 10, 10);
    mesh_drawable shank = mesh_primitive_cylinder(0.015f, { 0,0,-0.2f }, { 0,0,0 }, 10, 10);
    mesh_drawable knee = mesh_primitive_sphere(0.03f, { 0,0,0 }, 10, 10);
    mesh_drawable foot = mesh_drawable(mesh_load_file_obj("scenes/3D_graphics/01_modeling/assets/camel_foot.obj"));
    mesh_drawable eye = mesh_primitive_sphere(0.01f, { 0,0,0 }, 10, 10);
    articulation_point.uniform.color = { 1.0f, 0.75f, 0.2f };
    articulation_point.uniform.shading.specular = 0.01f;
    thigh.uniform.color = { 1.0f, 0.75f, 0.2f };
    thigh.uniform.shading.specular = 0.01f;
    shank.uniform.color = { 1.0f, 0.75f, 0.2f };
    shank.uniform.shading.specular = 0.01f;
    knee.uniform.color = { 1.0f, 0.75f, 0.2f };
    knee.uniform.shading.specular = 0.01f;
    foot.uniform.color = { 1.0f, 0.75f, 0.2f };
    foot.uniform.shading.specular = 0.01f;
    foot.uniform.transform.rotation = vcl::rotation_from_axis_angle_mat3({ 0,0,1 }, 3.14f) * R_camel;
    foot.uniform.transform.scaling = 2.0f;
    eye.uniform.color = { 0.1f,0.1f,0.1f };
    // assembling camel parts:
    camel.add(camel_trunk, "trunk");
    camel.add(camel_head, "head", "trunk", { 0, -0.17f, 0 });
    camel.add(eye, "right_eye", "head", { -0.016f, -0.195f, 0.36f });
    camel.add(eye, "left_eye", "head", { 0.016f, -0.195f, 0.36f });
    camel.add(articulation_point, "articulation_point_back_right", "trunk", { -0.08f, 0.1f, -0.1f });
    camel.add(articulation_point, "articulation_point_back_left", "trunk", { 0.08f, 0.1f, -0.1f });
    camel.add(articulation_point, "articulation_point_front_right", "trunk", { -0.08f, -0.1f, -0.1f });
    camel.add(articulation_point, "articulation_point_front_left", "trunk", { 0.08f, -0.1f, -0.1f });
    camel.add(thigh, "thigh_back_right", "articulation_point_back_right", { 0,0,0 });
    camel.add(thigh, "thigh_back_left", "articulation_point_back_left", { 0,0,0 });
    camel.add(thigh, "thigh_front_right", "articulation_point_front_right", { 0,0,0 });
    camel.add(thigh, "thigh_front_left", "articulation_point_front_left", { 0,0,0 });
    camel.add(knee, "knee_back_right", "thigh_back_right", { 0,0,-0.2f });
    camel.add(knee, "knee_back_left", "thigh_back_left", { 0,0,-0.2f });
    camel.add(knee, "knee_front_right", "thigh_front_right", { 0,0,-0.2f });
    camel.add(knee, "knee_front_left", "thigh_front_left", { 0,0,-0.2f });
    camel.add(shank, "shank_back_right", "knee_back_right", { 0,0,0 });
    camel.add(shank, "shank_back_left", "knee_back_left", { 0,0,0 });
    camel.add(shank, "shank_front_right", "knee_front_right", { 0,0,0 });
    camel.add(shank, "shank_front_left", "knee_front_left", { 0,0,0 });
    camel.add(foot, "foot_back_right", "shank_back_right", { 0,0,-0.19f });
    camel.add(foot, "foot_back_left", "shank_back_left", { 0,0,-0.19f });
    camel.add(foot, "foot_front_right", "shank_front_right", { 0,0,-0.19f });
    camel.add(foot, "foot_front_left", "shank_front_left", { 0,0,-0.19f });
    camel.set_shader_for_all_elements(shaders["mesh"]);
    // camel settings   

    // ****** tree ******************
    tree = create_tree();
    tree2 = create_tree();
    
    // Display Skybox
    skybox = mesh_drawable(mesh_skybox());
    skybox.uniform.shading.specular = 0;
    skybox.texture_id = create_texture_gpu(image_load_png("scenes/3D_graphics/01_modeling/sunset/sunset_hipshot.png"));
    //----------------
    tree.set_shader_for_all_elements(shaders["mesh"]);
    tree2.set_shader_for_all_elements(shaders["mesh"]);

    //timer.scale = 0.5f; // speed in which t varies
    timer.t_max = 10.0f; // t goes from 0 to 10 and restart in 0
}



/** This function is called at each frame of the animation loop.
    It is used to compute time-varying argument and perform data data drawing */
void scene_model::frame_draw(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& )
{
    timer.update();
    set_gui();

    gui_scene.wireframe = false;

    const float t = timer.t;

    glEnable( GL_POLYGON_OFFSET_FILL ); // avoids z-fighting when displaying wireframe
    
    // DISPLAYING SKYBOX -- must be rendered first in order to work
    glBindTexture(GL_TEXTURE_2D, skybox.texture_id);
    glDepthMask(GL_FALSE);// https://learnopengl.com/Advanced-OpenGL/Cubemaps
    skybox.uniform.transform.translation = scene.camera.camera_position() + vcl::vec3(-0.5f, -0.5f, -0.5f);
    draw(skybox, scene.camera, shaders["mesh"]);
    glBindTexture(GL_TEXTURE_2D, scene.texture_white);
    glDepthMask(GL_TRUE);
                                        
    // DISPLAYING THE TERRAIN
    // Before displaying a textured surface: bind the associated texture id
    glBindTexture(GL_TEXTURE_2D, terrain.texture_id);
    //---------------------------------
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);
    //----------------------------------
    // Display terrain
    glPolygonOffset( 1.0, 1.0 );
    draw(terrain, scene.camera, shaders["mesh"]);
    //------------------------------------
    // After the surface is displayed it is safe to set the texture id to a white image
    //  Avoids to use the previous texture for another object
    glBindTexture(GL_TEXTURE_2D, scene.texture_white);
    
    //DISPLAYING OTHER ELEMENTS OF THE SCENE:
    // Pyramid:
    draw(pyramid, scene.camera, shaders["mesh"]);
    // Camel:
    mat3 change_of_orientation = vcl::rotation_from_axis_angle_mat3({ 0,0,1 }, 3.14f / 12);
    mat3 small_inclination_matrix = vcl::rotation_from_axis_angle_mat3({ 1,0,0 }, 3.14f / 24);
    camel["trunk"].transform.rotation = small_inclination_matrix * change_of_orientation; // border of oasis
    camel["trunk"].transform.translation = { 4.4f, -2.5f, 0.48f - 0.4f * camel_position(t, timer.t_max)};
    camel["head"].transform.rotation = vcl::rotation_from_axis_angle_mat3({ 1,0,0 }, 5 * 3.14f / 12 * neck_position(t, timer.t_max));
    camel["shank_back_right"].transform.rotation = vcl::rotation_from_axis_angle_mat3({ 1,0,0 }, 3.14f * leg_position(t, timer.t_max));
    camel["shank_back_left"].transform.rotation = vcl::rotation_from_axis_angle_mat3({ 1,0,0 }, 3.14f * leg_position(t, timer.t_max));
    camel["shank_front_right"].transform.rotation = vcl::rotation_from_axis_angle_mat3({ 1,0,0 }, 3.14f * leg_position(t, timer.t_max));
    camel["shank_front_left"].transform.rotation = vcl::rotation_from_axis_angle_mat3({ 1,0,0 }, 3.14f * leg_position(t, timer.t_max));
    camel["thigh_back_right"].transform.rotation = vcl::rotation_from_axis_angle_mat3({ 1,0,0 }, -3.14f / 2 * leg_position(t, timer.t_max));
    camel["thigh_back_left"].transform.rotation = vcl::rotation_from_axis_angle_mat3({ 1,0,0 }, -3.14f / 2 * leg_position(t, timer.t_max));
    camel["thigh_front_right"].transform.rotation = vcl::rotation_from_axis_angle_mat3({ 1,0,0 }, -3.14f / 2 * leg_position(t, timer.t_max));
    camel["thigh_front_left"].transform.rotation = vcl::rotation_from_axis_angle_mat3({ 1,0,0 }, -3.14f / 2 * leg_position(t, timer.t_max));
    camel.update_local_to_global_coordinates();
    draw(camel, scene.camera);
    // Tree:
    tree["tronc"].transform.translation = { 4.6f, -2.6f, -0.2f };
    tree.update_local_to_global_coordinates();
    draw(tree, scene.camera);
    tree2["tronc"].transform.rotation = vcl::rotation_from_axis_angle_mat3({ 0,0,1 }, 3.14f / 2);
    tree2["tronc"].transform.translation = { 4.3f, -4.8f, -0.2f };
    tree2.update_local_to_global_coordinates();
    draw(tree2, scene.camera);
    
    //------------------------------------
    if (gui_scene.wireframe) { // wireframe if asked from the GUI
        glPolygonOffset(1.0, 1.0);
        draw(terrain, scene.camera, shaders["wireframe"]);
        draw(pyramid, scene.camera, shaders["wireframe"]);
        draw(tree, scene.camera, shaders["wireframe"]);
    }

    // (Transparent element at the end)
    // DISPLAYING OASIS:
    water.uniform.transform.translation = { 6.f, -4.f, -0.1f };
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(false);
    draw(water, scene.camera, shaders["mesh"]);
    glDepthMask(true);


}

mesh mesh_skybox()
{
    mesh skybox;
    // Compute normal of the quadrangle (orthogonal to the two basis vectors p00->p01 and p00->p10)
    vec3 p[] = { {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, //bottom
                 {0,1,1}, {1,1,1}, {1,0,1}, {0,0,1}, //top
                 {0,1,1}, {0,0,1}, {0,0,0}, {0,1,0}, //left
                 {0,0,1}, {1,0,1}, {1,0,0}, {0,0,0}, //front
                 {1,0,1}, {1,1,1}, {1,1,0}, {1,0,0}, //right
                 {1,1,1}, {0,1,1}, {0,1,0}, {1,1,0} }; //back
    
    vec3 n[6];
    for (size_t ku = 0; ku < 6; ++ku) {
        n[ku] = normalize(cross(normalize(p[4*ku+1] - p[4*ku]), normalize(p[4*ku+3] - p[4*ku])));
        skybox.position.push_back(p[4*ku]);
        skybox.position.push_back(p[4*ku + 1]);
        skybox.position.push_back(p[4*ku + 2]);
        skybox.position.push_back(p[4*ku + 3]);
        skybox.normal.push_back(n[ku]);
        skybox.normal.push_back(n[ku]);
        skybox.normal.push_back(n[ku]);
        skybox.normal.push_back(n[ku]);
    }

    skybox.texture_uv = { {1/4.0f,2/3.0f}, {1/2.0f,2/3.0f}, {1/2.0f,1}, {1/4.0f,1}, //bottom
                          {1/4.0f,0}, {1/2.0f,0}, {1/2.0f,1/3.0f}, {1/4.0f,1/3.0f}, //top
                          {0,1/3.0f}, {1 / 4.0f,1 / 3.0f}, {1 / 4.0f,2 / 3.0f}, {0,2/3.0f}, //left
                          {1/4.0f,1/3.0f}, {1/2.0f,1/3.0f}, {1/2.0f,2/3.0f}, {1/4.0f,2/3.0f}, //front
                          {1/2.0f,1/3.0f}, {3/4.0f,1/3.0f}, {3/4.0f,2/3.0f}, {1/2.0f,2/3.0f}, //right
                          {3/4.0f,1/3.0f}, {1,1/3.0f}, {1,2/3.0f}, {3/4.0f,2/3.0f} }; //back

    for(size_t ku = 0; ku<6;++ku) {
        const uint3 triangle_1 = { 4 * ku,4 * ku + 1,4 * ku + 2 };
        const uint3 triangle_2 = { 4 * ku,4 * ku + 2,4 * ku + 3 };
        skybox.connectivity.push_back(triangle_1);
        skybox.connectivity.push_back(triangle_2);           // Quadrangle made up of two triangles
    }

    return skybox;
}

void build_dune(vec2 center, float M, float L, float u, float v, float& z)
{
    if (abs(u - center[0]) <= L / 2 && abs(v - center[1]) <= L / 2)
    {
        u = u - center[0] + L / 2;
        v = v - center[1] + L / 2;
        float v_curve = 2 / L * (u - L / 2) * (u - L / 2) + L / 2;
        float z_curve = M - (4 / (L * L)) * M * (v_curve - L / 2) * (v_curve - L / 2);
        if (v - v_curve > 0) {
            z += z_curve * (v - L) * (v - L) / ((L - v_curve) * (L - v_curve));
        }
        else {
            z += -z_curve * v * (v - 2 * v_curve) / (v_curve * v_curve);
        }
    }
}

// Evaluate height of the terrain for any (u,v) \in [0,1]
float evaluate_terrain_z(float u, float v)
{
    float z = 0;

    // Changed direction // I have used a rotation of 90 degrees 
    build_dune({ 0.25f, 0.25f }, 2.0f, 0.5, v, u, z);
    build_dune({ 0.5f, 0.25f }, 1.0f, 0.5, v, u, z);
    build_dune({ 0.4f, 0.4f }, 0.5f, 0.8, v, u, z);
    build_dune({ 0.75f, 0.25f }, 0.8f, 0.5, v, u, z);

    // oasis - hole for the lake
    float d = norm(vec2(u, v) - vec2(0.8, 0.3))/0.08f;
    z -= 0.5 * exp(-d * d);

    return z;
}

// Evaluate 3D position of the terrain for any (u,v) \in [0,1]
vec3 evaluate_terrain(float u, float v)
{
    const float x = 20*(u-0.5f);
    const float y = 20*(v-0.5f);
    const float z = evaluate_terrain_z(u,v);

    return {x,y,z};
}

// Generate terrain mesh
mesh create_terrain()
{
    // Number of samples of the terrain is N x N
    const size_t N = 1000;

    mesh terrain; // temporary terrain storage (CPU only)
    terrain.position.resize(N * N);
    terrain.color.resize(N * N);
    terrain.texture_uv.resize(N * N);

    // Fill terrain geometry
    for (size_t ku = 0; ku < N; ++ku)
    {
        for (size_t kv = 0; kv < N; ++kv)
        {
            // Compute local parametric coordinates (u,v) \in [0,1]
            const float u = ku / (N - 1.0f);
            const float v = kv / (N - 1.0f);

            // get gui parameters
            const float scaling = 1.2f;
            const int octave = 4;
            const float persistency = 0.5f;
            
            // Evaluate Perlin noise
            //0.8, 0.3
            float noise;
            if ((u-0.8)* (u - 0.8)+(v-0.3)* (v - 0.3)<0.01)
                noise = 1;
            else
                noise = perlin(scaling * u, scaling * v, octave, persistency);

            // 3D vertex coordinates
            vec3 p = evaluate_terrain(u, v);
            const float x = p.x;
            const float y = p.y;
            const float z = p.z * noise;

            const float c = 0.3f + 0.7f * noise;

            // Compute coordinates
            terrain.position[kv + N * ku] = { x,y,z };
            terrain.color[kv + N * ku] = { c,c,c,1.0f };

            // Add texture
            int tiles_param = 10;
            terrain.texture_uv[kv + N * ku] = { (float)tiles_param * ku / N,(float)tiles_param * kv / N };
        }
    }


    // Generate triangle organization
    //  Parametric surface with uniform grid sampling: generate 2 triangles for each grid cell
    const unsigned int Ns = N;
    for (unsigned int ku = 0; ku < Ns - 1; ++ku)
    {
        for (unsigned int kv = 0; kv < Ns - 1; ++kv)
        {
            const unsigned int idx = kv + N * ku; // current vertex offset

            const uint3 triangle_1 = { idx, idx + 1 + Ns, idx + 1 };
            const uint3 triangle_2 = { idx, idx + Ns, idx + 1 + Ns };

            terrain.connectivity.push_back(triangle_1);
            terrain.connectivity.push_back(triangle_2);
        }
    }

    return terrain;
}

mesh create_pyramid(float side, float height, float z_offset, float rot) {
    const vec3 p0 = { (side / (float)sqrt(2)) * std::cos(2 * 3.14f * (0.0f/4 - rot)), (side / (float)sqrt(2)) * std::sin(2 * 3.14f * (0.0f/4 - rot)), z_offset };
    const vec3 p1 = { (side / (float)sqrt(2)) * std::cos(2 * 3.14f * (1.0f/4 - rot)), (side / (float)sqrt(2)) * std::sin(2 * 3.14f * (1.0f/4 - rot)), z_offset };
    const vec3 p2 = { (side / (float)sqrt(2)) * std::cos(2 * 3.14f * (2.0f/4 - rot)), (side / (float)sqrt(2)) * std::sin(2 * 3.14f * (2.0f/4 - rot)), z_offset };
    const vec3 p3 = { (side / (float)sqrt(2)) * std::cos(2 * 3.14f * (3.0f / 4 - rot)), (side / (float)sqrt(2)) * std::sin(2 * 3.14f * (3.0f / 4 - rot)), z_offset };
    const vec3 p4 = { 0,0,height + z_offset };

    mesh shape;
    shape.position = { p0, p1, p2, p3,
                       p0, p1, p4,
                       p1, p2, p4,
                       p2, p3, p4,
                       p3, p0, p4
                     };


    const vec3 n1 = normalize((p1-p0)*(p3-p0));
    const vec3 n2 = normalize((p1 - p0) * (p4 - p0));
    const vec3 n3 = normalize((p2 - p1) * (p4 - p1));
    const vec3 n4 = normalize((p3-p2)*(p4-p2));
    const vec3 n5 = normalize((p0-p3)*(p4-p3));

    shape.normal = { -n1, -n1, -n1, -n1,
                     -n2, -n2, -n2,
                     n3, n3, n3,
                     n4, n4, n4,
                     n5, n5, n5
                     };

    shape.connectivity = { {0,1,2}, {0,2,3},
                          {4,5,6},
                          {7,8,9},
                          {10,11,12},
                          {13,14,15}
                         };

    return shape;
}


mesh create_tronc(float radius, float height)
{
    mesh m;

    // Number of samples
    const size_t N = 20;

    // Number of divisions
    const size_t M = 7;

    // Geometry
    for (size_t k = 0; k < N; ++k)
    {
        const float u = k / float(N);
        const vec3 p = { radius * std::cos(2 * 3.14f * u), radius * std::sin(2 * 3.14f * u), 0.0f };
        m.position.push_back(p);
        //m.position.push_back(p + vec3(0, 0, height));
        for (size_t w = 1; w < M; ++w)
        {
            const vec3 p = { radius* (1-((float)0.5*w / (M - 1))) * std::cos(2 * 3.14f * u), radius * (1- ((float)0.5*w/(M-1))) * std::sin(2 * 3.14f * u), 0.0f };
            //std::cout << p.x << " " << p.y << std::endl;
            m.position.push_back(p +vec3(3*radius*std::sqrt((float)w/(M-1)), 0, w * height / (M - 1)));
        }
    }

    // Connectivity
    for (size_t k = 0; k < N; ++k)
    {
        for (size_t w = 0; w < M-1; ++w) 
        {
            const unsigned int u00 = M * k + w;
            const unsigned int u01 = (M * k + w + 1) % (M * N);
            const unsigned int u10 = (M * (k + 1)+w) % (M * N);
            const unsigned int u11 = (M * (k + 1) + w + 1) % (M * N);

            const uint3 t1 = { u00, u10, u11 };
            const uint3 t2 = { u00, u11, u01 };
            m.connectivity.push_back(t1);
            m.connectivity.push_back(t2);
        }
        
    }

    return m;
}

mesh create_foliage(float radius, float height, float z_offset)
{
    mesh m;

    // Cocotiers folliage
    // *************************** //

    const size_t N = 15; // Feuilles in each side

    // geometry
    const float d_2 = radius / float(2*N);
    for (size_t k = 0; k < N; ++k)
    {
        const vec3 p1 = { 2*k*d_2, 0.0f, -z_offset *k*k/float(N*N)};
        const vec3 p2 = { (2 * k + 1) * d_2, 5 * z_offset * ((k+0.5f) / float(N) - (k+0.5f) * (k+0.5f) / float(N * N)), -z_offset };
        const vec3 p3 = { (2 * k + 1) * d_2, -5 * z_offset * ((k + 0.5f) / float(N) - (k + 0.5f) * (k + 0.5f) / float(N * N)), -z_offset };
        m.position.push_back(p1 + vec3{ 0,0,height });
        m.position.push_back(p2 + vec3{ 0,0,height });
        m.position.push_back(p3 + vec3{ 0,0,height });
    }
    // apex
    m.position.push_back({ radius,0,height - z_offset });

    // connectivity
    const unsigned int Ns = N;
    for (unsigned int k = 0; k < Ns; ++k) {
        m.connectivity.push_back({ 3 * k , 3 * k + 1, 3 * k + 3 });
        m.connectivity.push_back({ 3 * k , 3 * k + 2, 3 * k + 3 });
    }

    return m;
}

hierarchy_mesh_drawable create_tree()
{
    //Tronc cocotiers:
    float tree_r_param = 0.1f;
    float tree_height = 2.0f;
    mesh_drawable tronc = mesh_drawable(create_tronc(tree_r_param, tree_height));
    tronc.uniform.color = { 0.5f, 0.25f, 0 };
    //Folliage cocotiers:
    mesh_drawable folliage = mesh_drawable(create_foliage(1.0, 0, 0.1));
    folliage.uniform.color = { 0, 0.25f, 0 };
    // détail : haut du tronc:
    mesh_drawable haut = mesh_drawable(mesh_primitive_sphere(tree_r_param / 2, { 0,0,0 }, 20, 20));
    haut.uniform.color = { 0.5f, 0.25f, 0 };
    haut.uniform.shading.specular = tronc.uniform.shading.specular;
    haut.uniform.shading.diffuse = tronc.uniform.shading.diffuse;
    // cocos
    mesh_drawable coco = mesh_drawable(mesh_primitive_sphere(3.0f * tree_r_param / 4, { 0,0,0 }, 20, 20));
    coco.uniform.color = { 111.0f / 255,138.0f / 255,37.0f / 255 };
    coco.uniform.shading.specular = 0.1f;
    coco.uniform.shading.diffuse = tronc.uniform.shading.diffuse;
    //----------------
    mat3 R1 = vcl::rotation_from_axis_angle_mat3({ 0,0,1 }, 2 * 3.14f / 3);
    mat3 R2 = vcl::rotation_from_axis_angle_mat3({ 0,0,1 }, 2 * 3.14f / 6);
    mat3 R3 = vcl::rotation_from_axis_angle_mat3({ 0,0,1 }, 2 * 3.14f / 12);
    mat3 Ry = vcl::rotation_from_axis_angle_mat3({ 0,1,0 }, 2 * 3.14f / 24);
    //----------------
    hierarchy_mesh_drawable tree;
    tree.add(tronc, "tronc");
    tree.add(haut, "haut", "tronc", { 3 * tree_r_param,0,tree_height });
    tree.add(coco, "coco1", "tronc", { 1.5f * tree_r_param,0,9 * tree_height / 10 });
    tree.add(coco, "coco2", "tronc", { 4.5f * tree_r_param,0,9 * tree_height / 10 });
    tree.add(coco, "coco3", "tronc", { 3.0f * tree_r_param,1.5f * tree_r_param,9 * tree_height / 10 });
    tree.add(coco, "coco4", "tronc", { 3.0f * tree_r_param,-1.5f * tree_r_param,9 * tree_height / 10 });
    tree.add(coco, "coco5", "tronc", { 4.5f * tree_r_param,0,0 });
    tree.add(coco, "coco6", "tronc", { tree_r_param,-1.5f * tree_r_param,tree_r_param });
    tree.add(folliage, "feuille_1", "tronc", { 3 * tree_r_param,0,tree_height });
    tree.add(folliage, "feuille_2", "tronc", { { 3 * tree_r_param,0,tree_height } , R1 });
    tree.add(folliage, "feuille_3", "tronc", { { 3 * tree_r_param,0,tree_height } , R1 * R1 });
    tree.add(folliage, "feuille_4", "tronc", { { 3 * tree_r_param,0,tree_height } , R2 * Ry });
    tree.add(folliage, "feuille_5", "tronc", { { 3 * tree_r_param,0,tree_height } , R1 * R2 * Ry });
    tree.add(folliage, "feuille_6", "tronc", { { 3 * tree_r_param,0,tree_height } , R1 * R1 * R2 * Ry });
    tree.add(folliage, "feuille_7", "tronc", { { 3 * tree_r_param,0,tree_height } , R3 * Ry * Ry });
    tree.add(folliage, "feuille_8", "tronc", { { 3 * tree_r_param,0,tree_height } , R3 * R3 * R3 * Ry * Ry });
    tree.add(folliage, "feuille_9", "tronc", { { 3 * tree_r_param,0,tree_height } , R1 * R3 * Ry * Ry });
    tree.add(folliage, "feuille_10", "tronc", { { 3 * tree_r_param,0,tree_height } , R1 * R3 * R3 * R3 * Ry * Ry });
    tree.add(folliage, "feuille_11", "tronc", { { 3 * tree_r_param,0,tree_height } , R1 * R1 * R3 * Ry * Ry });
    tree.add(folliage, "feuille_12", "tronc", { { 3 * tree_r_param,0,tree_height } , R1 * R1 * R3 * R3 * R3 * Ry * Ry });

    return tree;
}

float camel_position(float t, float& t_max)
{
    if (t < t_max / 6) // 0 -> t_max/6: camel stands up
        return 0;

    if (t_max / 6 <= t && t < t_max / 3) // t_max/6 -> 2*t_max/6 : camel ITSELF is bending down
        return std::sin(3.14f/2 * std::sin((t - t_max / 6) * 3 * 3.14f / t_max) );

    if (t_max / 3 <= t && t < t_max / 2) // 2*t_max/6 -> 3*t_max/6 : camel's HEAD is bending down to drink water
        return 1;

    if (t_max / 2 <= t && t < 2 * t_max / 3) // 3*t_max/6 -> 4*t_max/6 : camel fixes his head down to drink water
        return 1;

    if (2 * t_max / 3 <= t && t < 5 * t_max / 6) // 4*t_max/6 -> 5*t_max/6 : camel's HEAD is bending up
        return 1;

    if (5 * t_max / 6 <= t && t <= t_max) // t_max/6 -> 2*t_max/6 : camel ITSELF is bending up
        return std::sin(3.14f/2 * std::sin((t - 5 * t_max / 6) * (-3) * 3.14f / t_max + 3.14f / 2) );

    return 0; // au cas ou
}

float neck_position(float t, float &t_max)
{
    if (t < t_max / 6) // 0 -> t_max/6: camel stands up
        return 0;
    
    if (t_max / 6 <= t && t < t_max / 3) // t_max/6 -> 2*t_max/6 : camel ITSELF is bending down
        return 0;
    
    if (t_max / 3 <= t && t < t_max / 2) // 2*t_max/6 -> 3*t_max/6 : camel's HEAD is bending down to drink water
        return std::sin((t - t_max / 3) * 3 * 3.14f / t_max);

    if (t_max / 2 <= t && t < 2 * t_max / 3) // 3*t_max/6 -> 4*t_max/6 : camel fixes his head down to drink water
        return 1;

    if (2 * t_max / 3 <= t && t < 5 * t_max / 6) // 4*t_max/6 -> 5*t_max/6 : camel's HEAD is bending up
        return std::sin((t - 2 * t_max / 3) * (-3) * 3.14f / t_max + 3.14f / 2);

    if (5 * t_max / 6 <= t && t <= t_max) // t_max/6 -> 2*t_max/6 : camel ITSELF is bending up
        return 0;

    return 0; // au cas ou
}

float leg_position(float t, float& t_max)
{
    if (t < t_max / 6) // 0 -> t_max/6: camel stands up
        return 0;

    if (t_max / 6 <= t && t < t_max / 3) // t_max/6 -> 2*t_max/6 : camel ITSELF is bending down
        return std::sin((t - t_max / 6) * 3 * 3.14f / t_max);

    if (t_max / 3 <= t && t < t_max / 2) // 2*t_max/6 -> 3*t_max/6 : camel's HEAD is bending down to drink water
        return 1;

    if (t_max / 2 <= t && t < 2 * t_max / 3) // 3*t_max/6 -> 4*t_max/6 : camel fixes his head down to drink water
        return 1;

    if (2 * t_max / 3 <= t && t < 5 * t_max / 6) // 4*t_max/6 -> 5*t_max/6 : camel's HEAD is bending up
        return 1;

    if (5 * t_max / 6 <= t && t <= t_max) // t_max/6 -> 2*t_max/6 : camel ITSELF is bending up
        return std::sin((t - 5 * t_max / 6) * (-3) * 3.14f / t_max + 3.14f / 2);

    return 0; // au cas ou
}

void scene_model::set_gui()
{
    ImGui::Checkbox("Wireframe", &gui_scene.wireframe);

    ImGui::Spacing();
    ImGui::SliderFloat("Time", &timer.t, timer.t_min, timer.t_max);
    ImGui::SliderFloat("Time scale", &timer.scale, 0.1f, 3.0f);
}



#endif

