
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
    water = mesh_drawable(mesh_primitive_quad({ -1.f, -1.f, 0 }, { -1.f, 1.f, 0 }, { 1.f, 1.f, 0 }, { 1.f, -1.f, 0 }));
    water.uniform.color = { 0.6f, 0.6f, 0.9f };
    water.uniform.color_alpha = 0.5f; // transparency
    water.uniform.shading.specular = 0.8f;

    // Pyramid :
    pyramid = mesh_drawable(create_pyramid(6.0f, 4.0f, 0.0f, 0.25f));
    pyramid.uniform.color = { 241 / 255.0f,175 / 255.0f,0.0f };
    pyramid.uniform.shading.specular = 0.0f;
    pyramid.uniform.shading.ambiant = 0.5f; // MUDEI ISSO PQ TAVA MT ESCURO MAS PODEMOS DEIXAR OUTRO VALOR SE PREFERIR
    pyramid.uniform.transform.translation = { 5,5,0 };


    // ****** tree ******************
    //Tronc cocotiers:
    float tree_r_param = 0.1f;
    float tree_height = 2.0f;
    tronc = mesh_drawable(create_tronc(tree_r_param, tree_height));
    //Folliage cocotiers;
    folliage = mesh_drawable(create_foliage(1.0, tree_height, 0.1));
    //----------------
    vcl::mat3 R1 = vcl::rotation_from_axis_angle_mat3({ 0,0,1 }, 2 * 3.14f / 3);
    vcl::mat3 R2 = vcl::rotation_from_axis_angle_mat3({ 0,0,1 }, 2 * 3.14f / 6);
    vcl::mat3 Ry = vcl::rotation_from_axis_angle_mat3({ 0,1,0 }, 2 * 3.14f / 24);
    //----------------
    tree.add(tronc, "tronc");
    tree.add(folliage, "feuille_1", "tronc", { 3 * tree_r_param,0,0 });
    tree.add(folliage, "feuille_2", "tronc", { { 3 * tree_r_param,0,0 } , R1 });
    tree.add(folliage, "feuille_3", "tronc", { { 3 * tree_r_param,0,0 } , R1 * R1 });
    tree.add(folliage, "feuille_4", "tronc", { { 3 * tree_r_param,0,0 } , R2 * Ry });
    tree.add(folliage, "feuille_5", "tronc", { { 3 * tree_r_param,0,0 } , R1 * R2 * Ry });
    tree.add(folliage, "feuille_6", "tronc", { { 3 * tree_r_param,0,0 } , R1 * R1 * R2 * Ry });
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
    //draw(tronc, scene.camera, shaders["mesh"]);
    //draw(folliage, scene.camera, shaders["mesh"]);
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
    water.uniform.transform.translation = { 0, 0, -0.1f };
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(false);
    draw(water, scene.camera, shaders["mesh"]);
    glDepthMask(true);


}



// Evaluate height of the terrain for any (u,v) \in [0,1]
float evaluate_terrain_z(float u, float v)
{
    // medium dunes:
    int n_dunes = 4;
    vec2 center[4] = { {0.1f, 0.1f} , { 0.4f, 0.2f } , {0.2f, 0.6f}, {0.4f, 0.8f} };
    float M[4] = { 0.4f , 0.7f, 0.5f, 0.5f};
    float L = 0.2f;
    float z = 0;
    for (int i = 0; i < n_dunes; i++) {
        if (abs(u - center[i][0]) <= L/2 && abs(v - center[i][1]) <= L/2)
        {
            u = u - center[i][0] + L/2;
            v = v - center[i][1] + L/2;
            float v_curve = 2/L * (u - L/2) * (u - L/2) + L/2;
            float z_curve = M[i] - ( 4 / (L * L) ) * M[i] * (v_curve - L/2) * (v_curve - L/2);
            if (v - v_curve > 0) {
                z += z_curve * (v - L) * (v - L) / ((L - v_curve) * (L - v_curve));
            }
            else {
                z += -z_curve * v * (v - 2 * v_curve) / (v_curve * v_curve);
            }
        }
    }
    // small dunes:
    int n_dunes2 = 6;
    vec2 center2[6] = { {0.2f, 0.3f} , { 0.1f, 0.4f } , {0.4f, 0.5f}, {0.1f, 0.8f}, {0.15f, 0.95f}, {0.25f, 0.9f} };
    float M2[6] = { 0.2f , 0.5f, 0.3f, 0.1f, 0.4f, 0.4f };
    float L2 = 0.1f;
    for (int i = 0; i < n_dunes2; i++) {
        if (abs(u - center2[i][0]) <= L2 / 2 && abs(v - center2[i][1]) <= L2 / 2)
        {
            u = u - center2[i][0] + L2 / 2;
            v = v - center2[i][1] + L2 / 2;
            float v_curve = 2 / L2 * (u - L2 / 2) * (u - L2 / 2) + L2 / 2;
            float z_curve = M2[i] - (4 / (L2 * L2)) * M2[i] * (v_curve - L2 / 2) * (v_curve - L2 / 2);
            if (v - v_curve > 0) {
                z += z_curve * (v - L2) * (v - L2) / ((L2 - v_curve) * (L2 - v_curve));
            }
            else {
                z += -z_curve * v * (v - 2 * v_curve) / (v_curve * v_curve);
            }
        }
    }
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

