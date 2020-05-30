
#include "modeling.hpp"


#ifdef SCENE_3D_GRAPHICS

// Add vcl namespace within the current one - Allows to use function from vcl library without explicitely preceeding their name with vcl::
using namespace vcl;




float evaluate_terrain_z(float u, float v);
vec3 evaluate_terrain(float u, float v);
mesh create_terrain();
mesh create_pyramid(float radius, float height, float z_offset, float rot);
mesh create_tronc(float radius, float height);
mesh create_foliage(float radius, float height, float z_offset);


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
    scene.camera.apply_rotation(0,0,0,1.2f);

    // oasis small lake :
    water = mesh_drawable(mesh_primitive_quad({ -3.f, -3.f, 0 }, { -3.f, 3.f, 0 }, { 3.f, 3.f, 0 }, { 3.f, -3.f, 0 }));
    water.uniform.color = { 0.6f, 0.6f, 0.9f };
    water.uniform.color_alpha = 0.5f; // transparency
    water.uniform.shading.specular = 0.8f;

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
    pyramid.uniform.color = { 241 / 255.0f,175 / 255.0f,0.0f };
    pyramid.uniform.shading.specular = 0.0f;
    pyramid.uniform.shading.ambiant = 0.5f; // MUDEI ISSO PQ TAVA MT ESCURO MAS PODEMOS DEIXAR OUTRO VALOR SE PREFERIR
    

    // Camel
    // main part of camel:
    mesh_drawable camel_trunk = mesh_drawable(mesh_load_file_obj("scenes/3D_graphics/01_modeling/assets/camelo_trunk.obj"));
    mesh_drawable camel_head = mesh_drawable(mesh_load_file_obj("scenes/3D_graphics/01_modeling/assets/camelo_head.obj")); ////
    camel_trunk.uniform.color = { 1.0f, 0.75f, 0.2f };
    camel_trunk.uniform.shading.specular = 0.01f;
    camel_head.uniform.color = { 1.0f, 0.75f, 0.2f }; ////
    camel_head.uniform.shading.specular = 0.01f; ////
    vcl::mat3 R_camel = vcl::rotation_from_axis_angle_mat3({ 1,0,0 }, 3.14f / 2);
    camel_trunk.uniform.transform.rotation = R_camel;
    camel_head.uniform.transform.rotation = R_camel; ////
    camel_head.uniform.transform.translation = { 0, 0.17f, 0 }; ////
    // leg:
    mesh_drawable articulation_point = mesh_primitive_sphere(0.02f, { 0,0,0 }, 10, 10);
    mesh_drawable thigh = mesh_primitive_cylinder(0.025f, { 0,0,-0.2f }, { 0,0,0 }, 10, 10);
    mesh_drawable shank = mesh_primitive_cylinder(0.015f, { 0,0,-0.2f }, { 0,0,0 }, 10, 10);
    mesh_drawable knee = mesh_primitive_sphere(0.03f, { 0,0,0 }, 10, 10);
    mesh_drawable foot = mesh_drawable(mesh_load_file_obj("scenes/3D_graphics/01_modeling/assets/camel_foot.obj"));
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
    // assembling camel parts:
    camel.add(camel_trunk, "trunk");
    camel.add(camel_head, "head", "trunk", { 0, -0.17f, 0 }); ////
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
    //Tronc cocotiers:
    float tree_r_param = 0.1f;
    float tree_height = 2.0f;
    tronc = mesh_drawable(create_tronc(tree_r_param, tree_height));
    tronc.uniform.color = { 0.5f, 0.25f, 0 };
    //Folliage cocotiers:
    folliage = mesh_drawable(create_foliage(1.0, 0, 0.1));
    folliage.uniform.color = { 0, 0.25f, 0 };
    // détail : haut du tronc:
    mesh_drawable haut = mesh_drawable(mesh_primitive_sphere(tree_r_param / 2, { 0,0,0 }, 20, 20));
    haut.uniform.color = { 0.5f, 0.25f, 0 };
    haut.uniform.shading.specular = tronc.uniform.shading.specular;
    haut.uniform.shading.diffuse = tronc.uniform.shading.diffuse;
    // cocos
    mesh_drawable coco = mesh_drawable(mesh_primitive_sphere(3.0f*tree_r_param/4, { 0,0,0 }, 20, 20));
    coco.uniform.color = { 111.0f/255,138.0f/255,37.0f/255 };
    coco.uniform.shading.specular = 0.1f;
    coco.uniform.shading.diffuse = tronc.uniform.shading.diffuse;
    //----------------
    vcl::mat3 R1 = vcl::rotation_from_axis_angle_mat3({ 0,0,1 }, 2 * 3.14f / 3);
    vcl::mat3 R2 = vcl::rotation_from_axis_angle_mat3({ 0,0,1 }, 2 * 3.14f / 6);
    vcl::mat3 R3 = vcl::rotation_from_axis_angle_mat3({ 0,0,1 }, 2 * 3.14f / 12);
    vcl::mat3 Ry = vcl::rotation_from_axis_angle_mat3({ 0,1,0 }, 2 * 3.14f / 24);
    vcl::mat3 Ry2 = vcl::rotation_from_axis_angle_mat3({ 0,1,0 }, -2 * 3.14f / 24);
    //----------------
    tree.add(tronc, "tronc");
    tree.add(haut, "haut", "tronc", { 3 * tree_r_param,0,tree_height });
    tree.add(coco, "coco1", "tronc", { 1.5f * tree_r_param,0,9 * tree_height / 10 });
    tree.add(coco, "coco2", "tronc", { 4.5f * tree_r_param,0,9 * tree_height / 10 });
    tree.add(coco, "coco3", "tronc", { 3.0f * tree_r_param,1.5f*tree_r_param,9 * tree_height / 10 });
    tree.add(coco, "coco4", "tronc", { 3.0f * tree_r_param,-1.5f * tree_r_param,9 * tree_height / 10 });
    tree.add(coco, "coco5", "tronc", { 4.5f * tree_r_param,0,0 });
    tree.add(coco, "coco6", "tronc", { tree_r_param,-1.5f*tree_r_param,tree_r_param });
    tree.add(folliage, "feuille_1", "tronc", { 3 * tree_r_param,0,tree_height });
    tree.add(folliage, "feuille_2", "tronc", { { 3 * tree_r_param,0,tree_height } , R1 });
    tree.add(folliage, "feuille_3", "tronc", { { 3 * tree_r_param,0,tree_height } , R1 * R1 });
    tree.add(folliage, "feuille_4", "tronc", { { 3 * tree_r_param,0,tree_height } , R2 * Ry });
    tree.add(folliage, "feuille_5", "tronc", { { 3 * tree_r_param,0,tree_height } , R1 * R2 * Ry });
    tree.add(folliage, "feuille_6", "tronc", { { 3 * tree_r_param,0,tree_height } , R1 * R1 * R2 * Ry });
    tree.add(folliage, "feuille_7", "tronc", { { 3 * tree_r_param,0,tree_height } , R3 * Ry * Ry});
    tree.add(folliage, "feuille_8", "tronc", { { 3 * tree_r_param,0,tree_height } , R3 * R3 * R3 * Ry2 });
    tree.add(folliage, "feuille_9", "tronc", { { 3 * tree_r_param,0,tree_height } , R1 * R3 * Ry * Ry });
    tree.add(folliage, "feuille_10", "tronc", { { 3 * tree_r_param,0,tree_height } , R1 * R3 * R3 * R3 * Ry2 });
    tree.add(folliage, "feuille_11", "tronc", { { 3 * tree_r_param,0,tree_height } , R1 * R1 * R3 * Ry * Ry });
    tree.add(folliage, "feuille_12", "tronc", { { 3 * tree_r_param,0,tree_height } , R1 * R1 * R3 * R3 * R3 * Ry2 });
    //----------------
    tree.set_shader_for_all_elements(shaders["mesh"]);

}



/** This function is called at each frame of the animation loop.
    It is used to compute time-varying argument and perform data data drawing */
void scene_model::frame_draw(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& )
{
    set_gui();

    gui_scene.wireframe = false;

    glEnable( GL_POLYGON_OFFSET_FILL ); // avoids z-fighting when displaying wireframe
    //---------------------------------
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
    //Display Elements of the scene:
    draw(pyramid, scene.camera, shaders["mesh"]);
    camel["trunk"].transform.translation = { 4.0f, -2.5f, 0.48f };
    camel["head"].transform.rotation = vcl::rotation_from_axis_angle_mat3({ 1,0,0 }, 3.14f / 3); ////
    camel.update_local_to_global_coordinates();
    draw(camel, scene.camera);
    tree["tronc"].transform.translation = { 4.6f, -2.6f, -0.2f };
    tree.update_local_to_global_coordinates();
    draw(tree, scene.camera);   
    //------------------------------------
    // After the surface is displayed it is safe to set the texture id to a white image
    //  Avoids to use the previous texture for another object
    glBindTexture(GL_TEXTURE_2D, scene.texture_white);
    //------------------------------------
    if (gui_scene.wireframe) { // wireframe if asked from the GUI
        glPolygonOffset(1.0, 1.0);
        draw(terrain, scene.camera, shaders["wireframe"]);
        draw(pyramid, scene.camera, shaders["wireframe"]);
        draw(tronc, scene.camera, shaders["wireframe"]);
    }
    

    // Display oasis small lake
    water.uniform.transform.translation = { 6.f, -4.f, -0.1f };
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(false);
    draw(water, scene.camera, shaders["mesh"]);
    glDepthMask(true);


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

    build_dune({ 0.05f, 0.95f }, 0.20f, 0.07, u, v, z);
    build_dune({ 0.15f, 0.95f }, 0.22f, 0.08, u, v, z);
    build_dune({ 0.25f, 0.90f }, 0.24f, 0.07, u, v, z);
    build_dune({ 0.10f, 0.80f }, 0.25f, 0.10, u, v, z);
    build_dune({ 0.40f, 0.80f }, 0.22f, 0.10, u, v, z);
    build_dune({ 0.25f, 0.75f }, 0.26f, 0.10, u, v, z);
    build_dune({ 0.08f, 0.65f }, 0.22f, 0.08, u, v, z);
    build_dune({ 0.45f, 0.65f }, 0.25f, 0.07, u, v, z);
    build_dune({ 0.20f, 0.60f }, 0.25f, 0.10, u, v, z);
    build_dune({ 0.35f, 0.60f }, 0.27f, 0.08, u, v, z);
    build_dune({ 0.08f, 0.55f }, 0.22f, 0.07, u, v, z);
    build_dune({ 0.15f, 0.50f }, 0.25f, 0.08, u, v, z);
    build_dune({ 0.40f, 0.50f }, 0.25f, 0.10, u, v, z);
    build_dune({ 0.25f, 0.45f }, 0.29f, 0.07, u, v, z);
    build_dune({ 0.10f, 0.40f }, 0.25f, 0.10, u, v, z);
    build_dune({ 0.30f, 0.40f }, 0.23f, 0.07, u, v, z);
    build_dune({ 0.45f, 0.35f }, 0.27f, 0.08, u, v, z);
    build_dune({ 0.20f, 0.30f }, 0.22f, 0.07, u, v, z);
    build_dune({ 0.30f, 0.30f }, 0.23f, 0.07, u, v, z);
    build_dune({ 0.40f, 0.20f }, 0.28f, 0.10, u, v, z);
    build_dune({ 0.25f, 0.15f }, 0.24f, 0.07, u, v, z);
    build_dune({ 0.10f, 0.10f }, 0.23f, 0.10, u, v, z);
    build_dune({ 0.30f, 0.05f }, 0.27f, 0.08, u, v, z);
    
    // ONE SINGLE BIG DUNE:
    /*
    const vec2 u0 = { 0.5f, 0.5f };
    const float h0 = 2.0f;
    const float M = 3;
    const float v_curve = 2 * (u - 0.5) * (u - 0.5) + 0.5;
    float z = 0;
    float z_curve = 0;

    if (v - v_curve > 0) {
        z_curve = M - 4 * M * (v_curve - 0.5) * (v_curve - 0.5);
        z = z_curve * (v - 1) * (v - 1) / ((1 - v_curve) * (1 - v_curve));
    }
    else {
        z_curve = M - 4 * M * (v_curve - 0.5) * (v_curve - 0.5);
        z = -z_curve * v * (v - 2 * v_curve) / (v_curve * v_curve);
    }*/

    // smooth irregularitues under the dunes :
    float d1 = norm(vec2(u, v) - vec2(0.15, 0.3)) / 0.2f;
    float d2 = norm(vec2(u, v) - vec2(0.1, 0.8)) / 0.1f;
    float d3 = norm(vec2(u, v) - vec2(0.3, 0.5)) / 0.2f;
    z += 0.5 * exp(-d1 * d1) + 0.3 * exp(-d2 * d2) + 0.6 * exp(-d3 * d3);

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
    terrain.texture_uv.resize(N * N);

    // Fill terrain geometry
    for (size_t ku = 0; ku < N; ++ku)
    {
        for (size_t kv = 0; kv < N; ++kv)
        {
            // Compute local parametric coordinates (u,v) \in [0,1]
            const float u = ku / (N - 1.0f);
            const float v = kv / (N - 1.0f);

            // Compute coordinates
            terrain.position[kv + N * ku] = evaluate_terrain(u, v);

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

/* mesh create_pyramid(float side, float height, float z_offset, float rot)
{
    mesh m;

    // Pyramidal structure
    // *************************** //

    const size_t N = 4;

    // geometry
    for (size_t k = 0; k < N; ++k)
    {
        const float u = k / float(N);
        const vec3 p = { (side/(float)sqrt(2)) * std::cos(2 * 3.14f * (u-rot)), (side / (float)sqrt(2)) * std::sin(2 * 3.14f * (u-rot)), 0.0f };
        m.position.push_back(p + vec3{ 0,0,z_offset });
    }
    // apex
    m.position.push_back({ 0,0,height + z_offset });

    // connectivity
    const unsigned int Ns = N;
    for (unsigned int k = 0; k < Ns; ++k) {
        m.connectivity.push_back({ k , (k + 1) % Ns, Ns });
    }

    // close the bottom of the cone
    // *************************** //

    // Geometry
    for (size_t k = 0; k < N; ++k)
    {
        const float u = k / float(N);
        const vec3 p = { (side / (float)sqrt(2)) * std::cos(2 * 3.14f * (u-rot)), (side / (float)sqrt(2)) * std::sin(2 * 3.14f * (u-rot)), 0.0f };
        m.position.push_back(p + vec3{ 0,0,z_offset });
    }
    // central position
    m.position.push_back({ 0,0,z_offset });

    // connectivity
    for (unsigned int k = 0; k < Ns; ++k)
        m.connectivity.push_back({ k + Ns + 1, (k + 1) % Ns + Ns + 1, 2 * Ns + 1 });

    return m;
}*/

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

void scene_model::set_gui()
{
    ImGui::Checkbox("Wireframe", &gui_scene.wireframe);
}



#endif

