/*
 * Graphics Engine
 *
 * File:     render_pgm.c
 * Author:   James Bray
 * Repo:     https://github.com/James-Bray19/ASCII-Graphics-Engine
 *
 * A simple 3D graphics engine that reads in one or more obj files 
 * and renders them as a pgm image.
 * 
 * LML repo:        https://github.com/jamesbray03/Lightweight-Matrix-Library
 * 4x4 matrices:    http://www.codinglabs.net/article_world_view_projection_matrix.aspx
 * Bresenham's:     https://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm
 * Scanline:        https://en.wikipedia.org/wiki/Scanline_rendering
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "lml.h"

#define PI 3.14159265358979323846

// --------------- Mesh Structure ---------------

typedef struct {
    Matrix *verts; // matrix where each row is a vertex
    Matrix *norms; // matrix where each row is a normal
    Matrix *faces; // matrix where each row is a face
    Matrix *face_norms; // matrix where each row is a face normal
} Mesh;

// --------------- Function Prototypes ---------------

/// @brief copies the contents of specified obj file into mesh struct
/// @param m destination mesh
/// @param f file to read from
/// @note recommended to use a normalised mesh, otherwise scaling may be required
/// @note ignores textures and materials for now
void load_obj_file(Mesh *m, FILE *f);

/// @brief applies a transformation matrix to a mesh
/// @param m mesh to transform
/// @param trans transformation matrix
void transform(Mesh *m, Matrix *trans);

/// @brief returns a perspective projection matrix
/// @param width screen width
/// @param height screen height
/// @param fov field of view
/// @param z_near near clipping plane
/// @param z_far far clipping plane
Matrix *get_projection(int width, int height, int fov, int z_near, int z_far);

/// @brief returns a translation matrix
Matrix *get_translation(double x, double y, double z);

/// @brief returns a rotation matrix
Matrix *get_rotation(double x, double y, double z);

/// @brief returns a scale matrix
Matrix *get_scale(double x, double y, double z);

/// @brief culls backfaces from a mesh
/// @param m mesh to cull
/// @note apply before projection
void cull_backfaces(Mesh *m);

/// @brief draws a line on the screen at the specified coordinates with the specified intensity
void draw_line(Matrix *screen, int x1, int y1, int x2, int y2, double intensity);

/// @brief draws a triangle on the screen at the specified coordinates with the specified intensity
void draw_triangle(Matrix *screen, int x1, int y1, int x2, int y2, int x3, int y3, double line_colour, double fill_colour);

/// @brief renders the mesh to the screen
/// @param m mesh to render
/// @param width screen width
/// @param height screen height
/// @param perspective 0 for orthographic, 1 for perspective
void render(Mesh *m, int width, int height, int perspective, char *filename);

// --------------- Main Loop ---------------

void main() {
    char running = 'y';
    while (running == 'y') {
        printf("\e[1;1H\e[2J");
        printf("\n--------------- SCENE SETUP ---------------\n");



        printf("\nExport Settings\n");
        char pgm_name[256]; int width, height;

        do { printf("\trender name: "); scanf("%s", pgm_name); } while (strlen(pgm_name) == 0);
        do { printf("\tscreen width: "); scanf("%d", &width); } while (width <= 0);
        do { printf("\tscreen height: "); scanf("%d", &height); } while (height <= 0);

        // get the export file path
        char pgm_path[356];
        strcpy(pgm_path, "pgm/"); strcat(pgm_path, pgm_name); strcat(pgm_path, ".pgm");



        printf("\nCamera Settings\n");
        int perspective; double z_near, z_far, fov;

        char default_values = 'n';
        do { printf("    use default camera? (y/n): "); scanf(" %c", &default_values); } 
        while (default_values != 'y' && default_values != 'n');

        if (default_values == 'y') {
            perspective = 1; z_near = 0.1; z_far = 10; fov = 60;
        } else {
            do { printf("\tperspective: "); scanf("%d", &perspective); } while (perspective < 0 || perspective > 1);
            do { printf("\tnear clipping plane: "); scanf("%lf", &z_near); } while (z_near <= 0);
            do { printf("\tfar clipping plane: "); scanf("%lf", &z_far); } while (z_far <= 0);
            do { printf("\tfield of view: "); scanf("%lf", &fov); } while (fov <= 0 || fov >= 180);
        }
        Matrix *projection = get_projection(width, height, fov, z_near, z_far);



        printf("\nObject Settings");
        Mesh **queue = malloc(sizeof(Mesh *)); int num_objects = 0;

        char add_another = 'y';
        while (add_another == 'y') {
            num_objects++;

            // create mesh object
            queue = realloc(queue, num_objects * sizeof(Mesh *));
            queue[num_objects - 1] = malloc(sizeof(Mesh));

            printf("\n    Object %d:\n", num_objects);

            // get object file
            FILE *obj_file; char obj_name[256], obj_path[256];
            do {
                printf("\tobject: "); scanf("%s", obj_name);
                strcpy(obj_path, "obj/"); strcat(obj_path, obj_name); strcat(obj_path, ".obj");
                obj_file = fopen(obj_path, "r");
            } while (obj_file == NULL);
            load_obj_file(queue[num_objects - 1], obj_file); fclose(obj_file);

            // get transformation values
            char default_values = 'n';
            do { printf("\tuse default position? (y/n): "); scanf(" %c", &default_values); } 
            while (default_values != 'y' && default_values != 'n');

            double sx, sy, sz, rx, ry, rz, tx, ty, tz;
            if (default_values == 'y') {
                sx = 1; sy = 1; sz = 1;
                rx = 0; ry = 0; rz = 0;
                tx = 0; ty = 0; tz = 3;
            } else {
                printf("\t    scale: "); scanf("%lf %lf %lf", &sx, &sy, &sz);
                printf("\t    rotation: "); scanf("%lf %lf %lf", &rx, &ry, &rz);
                printf("\t    translation: "); scanf("%lf %lf %lf", &tx, &ty, &tz);
            }

            // position object in 3D space
            transform(queue[num_objects - 1], get_scale(sx, sy, sz));
            transform(queue[num_objects - 1], get_rotation(rx, ry, rz));
            transform(queue[num_objects - 1], get_translation(tx, ty, tz));
            cull_backfaces(queue[num_objects - 1]);

            // move to 2D space
            transform(queue[num_objects - 1], projection);

            do {printf("\tadd another object? (y/n): "); scanf(" %c", &add_another); }
            while (add_another != 'y' && add_another != 'n');
        }



        printf("\n--------------- FINISHED SETUP ---------------\n");
        printf("\nPress any key to render..."); getchar(); getchar();



        // create the pgm file
        FILE *pgm_file = fopen(pgm_path, "w");
        fprintf(pgm_file, "P2\n%d %d\n255\n", width, height);
        fclose(pgm_file);

        // render objects into the file
        for (int i = 0; i < num_objects; i++) { render(queue[i], width, height, perspective, pgm_path); }
        printf("\nScene rendered to %s\n", pgm_path);
        


        // clean up
        for (int i = 0; i < num_objects; i++) {
            release(queue[i]->verts); release(queue[i]->norms);
            release(queue[i]->faces); release(queue[i]->face_norms);
            free(queue[i]);
        }
    
        // ask to run again
        printf("\nRender another object? (y/n): ");
        scanf(" %c", &running);
    }
}

// --------------- Function Definitions ---------------

void load_obj_file(Mesh *m, FILE *f) {
    
    // count number of vertices and faces
    char line[256];
    int num_verts = 0, num_norms = 0, num_faces = 0;
    while (fgets(line, sizeof(line), f)) {

        // get the data type of the line
        if (line[0] == 'v') {
            if (line[1] == ' ') { num_verts++; }     // if data is a vertex
            else if (line[1] == 'n'){ num_norms++; } // if data is a normal
        } else if (line[0] == 'f') { num_faces++; }  // if data is a face
    }

    // create matrices to store vertices and faces
    m->verts = zeros(num_verts, 4);
    m->norms = zeros(num_norms, 3);
    m->faces = zeros(num_faces, 3);
    m->face_norms = zeros(num_faces, 3);

    // go back to the start
    rewind(f);

    // read in vertices and faces
    int vert_idx = 0, norm_idx = 0, face_idx = 0;
    while (fgets(line, sizeof(line), f)) {

        // if line is a vector
        if (line[0] == 'v') {

            // if line is a vertex
            if (line[1] == ' ') {
                sscanf(line, "v %lf %lf %lf", &m->verts->data[vert_idx][0], 
                                              &m->verts->data[vert_idx][1], 
                                              &m->verts->data[vert_idx][2]);
                m->verts->data[vert_idx][3] = 1; // sneaky homogenous coordinate
                vert_idx++;
            
            // if line is a normal
            } else if (line[1] == 'n') {
                sscanf(line, "vn %lf %lf %lf", &m->norms->data[norm_idx][0], 
                                               &m->norms->data[norm_idx][1], 
                                               &m->norms->data[norm_idx][2]);
                norm_idx++;
            }
        } 
        
        // if line is a face (group of indices)
        else if (line[0] == 'f') {

            // read in vertices and normals, ignore textures for now
            sscanf(line, "f %lf//%lf %lf//%lf %lf//%lf", 
            &m->faces->data[face_idx][0], &m->face_norms->data[face_idx][0],
            &m->faces->data[face_idx][1], &m->face_norms->data[face_idx][1],
            &m->faces->data[face_idx][2], &m->face_norms->data[face_idx][2]);
            face_idx++;
        }
    }
    fclose(f);
}

void transform(Mesh *m, Matrix *trans) {

    // placeholders for iteration
    Matrix *vert, *norm, *norm_trans;

    // transform all the vertices
    for (int i = 0; i < m->verts->rows; i++) {

        // get vertex
        vert = get_row(m->verts, i);
        vert = transpose(vert);
        vert = multiply(trans, vert);

        // divide by w to convert to 3D perpective
        if (vert->data[3][0] != 0) { vert = scalar_multiply(vert, 1 / vert->data[3][0]); }

        // then set the vertex
        vert = transpose(vert);
        set_row(m->verts, i, vert);
    }

    // transform all the normals
    for (int i = 0; i < m->norms->rows; i++) {

        // get normal
        norm = get_row(m->norms, i);
        norm = transpose(norm);

        // exclude the translation component of the 4x4
        norm_trans = get_submatrix(trans, 0, 0, 3, 3);
        
        // transform and normalise
        norm = multiply(norm_trans, norm);
        double mag = sqrt(norm->data[0][0] * norm->data[0][0] + 
                          norm->data[1][0] * norm->data[1][0] + 
                          norm->data[2][0] * norm->data[2][0]);
        norm = scalar_multiply(norm, 1 / mag);

        // and then set the normal
        norm = transpose(norm);
        set_row(m->norms, i, norm);
    }

    // release memory
    release(vert); release(norm); release(norm_trans);
}

Matrix *get_projection(int width, int height, int fov, int z_near, int z_far) {

    // pre-calculate common factors
    double tan_fov = tan(((fov/2) * PI) / 180.0);
    double aspect_ratio = (double)width / (double)height;

    // create projection matrix
    Matrix *projection = zeros(4, 4);
    projection->data[0][0] = 1 / (tan_fov * aspect_ratio);
    projection->data[1][1] = 1 / tan_fov;
    projection->data[2][2] = - (z_far + z_near) / (z_far - z_near);
    projection->data[2][3] = -1;
    projection->data[3][2] = - (2 * z_far * z_near) / (z_far - z_near);
    return projection;
}

Matrix *get_translation(double x, double y, double z) {

    // create translation matrix
    Matrix *translation = identity(4);
    translation->data[0][3] = x;
    translation->data[1][3] = y;
    translation->data[2][3] = z;
    return translation;
}

Matrix *get_rotation(double x, double y, double z) {

    // convert angles to radians
    double x_rad = x * PI / 180.0;
    double y_rad = y * PI / 180.0;
    double z_rad = z * PI / 180.0;

    // create rotation matrix
    Matrix *rotation = identity(4);
    rotation->data[0][0] = cos(z_rad) * cos(y_rad);
    rotation->data[0][1] = cos(z_rad) * sin(y_rad) * sin(x_rad) - sin(z_rad) * cos(x_rad);
    rotation->data[0][2] = cos(z_rad) * sin(y_rad) * cos(x_rad) + sin(z_rad) * sin(x_rad);
    rotation->data[1][0] = sin(z_rad) * cos(y_rad);
    rotation->data[1][1] = sin(z_rad) * sin(y_rad) * sin(x_rad) + cos(z_rad) * cos(x_rad);
    rotation->data[1][2] = sin(z_rad) * sin(y_rad) * cos(x_rad) - cos(z_rad) * sin(x_rad);
    rotation->data[2][0] = -sin(y_rad);
    rotation->data[2][1] = cos(y_rad) * sin(x_rad);
    rotation->data[2][2] = cos(y_rad) * cos(x_rad);
    return rotation;
}

Matrix *get_scale(double x, double y, double z) {

    // create scale matrix
    Matrix *scale = identity(4);
    scale->data[0][0] = x;
    scale->data[1][1] = y;
    scale->data[2][2] = z;
    return scale;
}

void cull_backfaces(Mesh *m) {

    // iterate backward so removing faces doesn't affect the loop
    for (int i = m->faces->rows - 1; i >= 0; i--) {

        // get the 'position' and normal of the face
        int vert_idx = m->faces->data[i][0] - 1; 
        int norm_idx = m->face_norms->data[i][0] - 1;
        Matrix *face_norm = get_row(m->norms, norm_idx);
        Matrix *face_vert = get_row(m->verts, m->faces->data[i][0] - 1);

        // determine if face is visible (if camera is moving use relative position here)
        double dot = face_norm->data[0][0] * (face_vert->data[0][0] - 0) +
                     face_norm->data[0][1] * (face_vert->data[0][1] - 0) +
                     face_norm->data[0][2] * (face_vert->data[0][2] - 0);
        
        // if unseen, cull the face
        if (dot > 0) {
            remove_row(m->faces, i);
            remove_row(m->face_norms, i);
        }
    }
}

double get_light_intensity(Mesh *m, int face_idx) {

    // get the normal of the face
    int norm_idx = m->face_norms->data[face_idx][0] - 1;
    Matrix *face_norm = get_row(m->norms, norm_idx);

    // get the light intensity
    double dot = face_norm->data[0][0] * 0.408248 +
                 face_norm->data[0][1] * 0.816497 +
                 face_norm->data[0][2] * 0.408248;

    return fmax(0.05, dot);
}

void draw_line(Matrix *screen, int x1, int y1, int x2, int y2, double intensity) {

    // Bresenham's line algorithm (see header)

    // calculate deltas, which direction to step and error
    int dx = abs(x2 - x1), dy = abs(y2 - y1);
    int sx = x1 < x2 ? 1 : -1, sy = y1 < y2 ? 1 : -1;
    int err = dx - dy, e2, done = 0;

    // loop until we reach the end
    while (!done) {

        // set pixel
        if (x1 >= 0 && x1 < screen->cols && y1 >= 0 && y1 < screen->rows) { 
            screen->data[y1][x1] = intensity; }
        
        // check if we've reached the end
        if (x1 == x2 && y1 == y2) { done = 1; } 

        // if not, calculate next pixel
        else {
            e2 = 2 * err;
            if (e2 > -dy) { err -= dy; x1 += sx; }
            if (e2 < dx) { err += dx; y1 += sy; }
        }
    }
}

void draw_triangle(Matrix *screen, int x1, int y1, int x2, int y2, int x3, int y3, double line_colour, double fill_colour) {
    
    // Scanline algorithm (see header)

    // sort vertices by vertical position
    int temp;
    if (y1 > y2) { temp = x1; x1 = x2; x2 = temp; temp = y1; y1 = y2; y2 = temp; }
    if (y1 > y3) { temp = x1; x1 = x3; x3 = temp; temp = y1; y1 = y3; y3 = temp; }
    if (y2 > y3) { temp = x2; x2 = x3; x3 = temp; temp = y2; y2 = y3; y3 = temp; }

    // calculate slopes
    double m1 = (y2 - y1 != 0) ? (double)(x2 - x1) / (double)(y2 - y1) : 0;
    double m2 = (y3 - y1 != 0) ? (double)(x3 - x1) / (double)(y3 - y1) : 0;
    double m3 = (y3 - y2 != 0) ? (double)(x3 - x2) / (double)(y3 - y2) : 0;

    // draw top half of triangle using scanline
    for (int y = y1; y <= y2; y++) {
        int x_start = x1 + m1 * (y - y1);
        int x_end = x1 + m2 * (y - y1);
        draw_line(screen, x_start, y, x_end, y, fill_colour);
    }

    // draw bottom half of triangle using scanline
    for (int y = y2; y <= y3; y++) {
        int x_start = x2 + m3 * (y - y2);
        int x_end = x1 + m2 * (y - y1);
        draw_line(screen, x_start, y, x_end, y, fill_colour);
    }

    // draw_line(screen, x1, y1, x2, y2, line_colour);
    // draw_line(screen, x2, y2, x3, y3, line_colour);
    // draw_line(screen, x3, y3, x1, y1, line_colour);
}

void buffer_to_pgm(Matrix *screen, int width, int height, char *filepath) {

    // read in previous image data and store in a buffer
    FILE *pgm_file = fopen(filepath, "r");
    if (pgm_file == NULL) { printf("Error opening file!\n"); return; }

    Matrix *prev_screen = zeros(height, width);

    // if file is empty leave prev screen as zeros
    if (fscanf(pgm_file, "P2\n%d %d\n255\n", &width, &height) == EOF) { fclose(pgm_file); }
    // otherwise fetch the previous screen data
    else {
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                fscanf(pgm_file, "%lf", &prev_screen->data[i][j]);
            }
        }
        fclose(pgm_file);
    }

    // rewrite image with current screen buffer (add object like stickers on top of previous image)
    pgm_file = fopen(filepath, "w");
    if (pgm_file == NULL) { printf("Error opening file!\n"); return; }
    fprintf(pgm_file, "P2\n%d %d\n255\n", width, height);
    
    Matrix *curr_screen = scalar_multiply(screen, 255);

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (screen->data[i][j] != 0) { 
                fprintf(pgm_file, "%d ", (int)curr_screen->data[i][j]); 
            } 
            else { 
                fprintf(pgm_file, "%d ", (int)prev_screen->data[i][j]); 
            }
        }
        fprintf(pgm_file, "\n");
    }
    fclose(pgm_file);
}

void render(Mesh *m, int width, int height, int perspective, char *filepath) {

    // create screen buffer
    Matrix *screen = zeros(height, width);
  
    // convert each vertex to screen space
    for (int i = 0; i < m->verts->rows; i++) {

        // apply perspective if selected
        if (perspective && m->verts->data[i][2] != 0) {
            m->verts->data[i][0] = m->verts->data[i][0] / m->verts->data[i][2];
            m->verts->data[i][1] = m->verts->data[i][1] / m->verts->data[i][2];
        }
        
        // normalise to screen space
        m->verts->data[i][0] = round((m->verts->data[i][0] + 1) * width / 2);
        m->verts->data[i][1] = round((m->verts->data[i][1] + 1) * height / 2);
    }

    // draw triangles and edges
    for (int i = 0; i < m->faces->rows; i++) {

        // get face vertices
        int v1 = m->faces->data[i][0] - 1;
        int v2 = m->faces->data[i][1] - 1;
        int v3 = m->faces->data[i][2] - 1;

        // get vertex coordinates
        int x1 = m->verts->data[v1][0]; int y1 = m->verts->data[v1][1];
        int x2 = m->verts->data[v2][0]; int y2 = m->verts->data[v2][1];
        int x3 = m->verts->data[v3][0]; int y3 = m->verts->data[v3][1];

        // fill triangle and draw edges (in rendering order)
        draw_triangle(screen, x1, y1, x2, y2, x3, y3, 1, get_light_intensity(m, i));
    }

    // write the buffer to a pgm file
    buffer_to_pgm(screen, width, height, filepath);
}