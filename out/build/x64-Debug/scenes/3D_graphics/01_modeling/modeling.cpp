
#include "modeling.hpp"


#ifdef SCENE_3D_GRAPHICS

// Add vcl namespace within the current one - Allows to use function from vcl library without explicitely preceeding their name with vcl::
using namespace vcl;




float evaluate_terrain_z(float u, float v);
vec3 evaluate_terrain(float u, float v);
mesh create_terrain();


/** This function is called before the beginning of the animation loop
    It is used to initialize all part-specific data */
void scene_model::setup_data(std::map<std::string,GLuint>& , scene_structure& scene, gui_structure& )
{
    // Create visual terrain surface
    terrain = create_terrain();
    terrain.uniform.color = {0.6f,0.85f,0.5f};
    terrain.uniform.shading.specular = 0.0f; // non-specular terrain material

    // Load a texture image on GPU and stores its ID
    texture_id = create_texture_gpu(image_load_png("scenes/3D_graphics/01_modeling/assets/sand.jpg"));

    // Setup initial camera mode and position
    scene.camera.camera_type = camera_control_spherical_coordinates;
    scene.camera.scale = 10.0f;
    scene.camera.apply_rotation(0,0,0,1.2f);

}



/** This function is called at each frame of the animation loop.
    It is used to compute time-varying argument and perform data data drawing */
void scene_model::frame_draw(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& )
{
    set_gui();

    glEnable( GL_POLYGON_OFFSET_FILL ); // avoids z-fighting when displaying wireframe

    // Before displaying a textured surface: bind the associated texture id
    glBindTexture(GL_TEXTURE_2D, texture_id);

    // Display terrain
    glPolygonOffset( 1.0, 1.0 );
    draw(terrain, scene.camera, shaders["mesh"]);

   // After the surface is displayed it is safe to set the texture id to a white image
   //  Avoids to use the previous texture for another object
    glBindTexture(GL_TEXTURE_2D, scene.texture_white);

    if( gui_scene.wireframe ){ // wireframe if asked from the GUI
        glPolygonOffset( 1.0, 1.0 );
        draw(terrain, scene.camera, shaders["wireframe"]);
    }
}



// Evaluate height of the terrain for any (u,v) \in [0,1]
float evaluate_terrain_z(float u, float v)
{
    const vec2 u0 = {0.5f, 0.5f};
    const float h0 = 2.0f;
    const float M = 3;
    const float v_curve = 2 * (u - 0.5) * (u - 0.5) + 0.5;
    float z = 0;
    float z_curve = 0;

    if (v - v_curve > 0) {
        z_curve = M - 4 * M * (v_curve - 0.5) * (v_curve - 0.5);
        z = z_curve * (v - 1) * (v - 1) / ((1 - v_curve) * (1 - v_curve));
    } else {
        z_curve = M - 4 * M * (v_curve - 0.5) * (v_curve - 0.5);
        z = -z_curve * v * (v - 2 * v_curve)/(v_curve*v_curve);
    }
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
    const size_t N = 100;

    mesh terrain; // temporary terrain storage (CPU only)
    terrain.position.resize(N*N);

    // Fill terrain geometry
    for(size_t ku=0; ku<N; ++ku)
    {
        for(size_t kv=0; kv<N; ++kv)
        {
            // Compute local parametric coordinates (u,v) \in [0,1]
            const float u = ku/(N-1.0f);
            const float v = kv/(N-1.0f);

            // Compute coordinates
            terrain.position[kv+N*ku] = evaluate_terrain(u,v);
        }
    }


    // Generate triangle organization
    //  Parametric surface with uniform grid sampling: generate 2 triangles for each grid cell
    const unsigned int Ns = N;
    for(unsigned int ku=0; ku<Ns-1; ++ku)
    {
        for(unsigned int kv=0; kv<Ns-1; ++kv)
        {
            const unsigned int idx = kv + N*ku; // current vertex offset

            const uint3 triangle_1 = {idx, idx+1+Ns, idx+1};
            const uint3 triangle_2 = {idx, idx+Ns, idx+1+Ns};

            terrain.connectivity.push_back(triangle_1);
            terrain.connectivity.push_back(triangle_2);
        }
    }

    return terrain;
}

void scene_model::set_gui()
{
    ImGui::Checkbox("Wireframe", &gui_scene.wireframe);
}



#endif

