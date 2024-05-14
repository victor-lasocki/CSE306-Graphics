#define _CRT_SECURE_NO_WARn_iNGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <cmath>
#include <iostream>
#include <chrono>
#include <list>
#include <random>

double square(double x) { return x * x;}

static std::default_random_engine engine(10);
static std::uniform_real_distribution<double> uniform ( 0 , 1 );

class Vector {
	public:
		explicit Vector(double x = 0, double y = 0, double z = 0) {
			data[0] = x;
			data[1] = y;
			data[2] = z;
		}
		double norm2() const {
			return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
		}
		double norm() const {
			return sqrt(norm2());
		}
		void normalize() {
			double n = norm();
			data[0] /= n;
			data[1] /= n;
			data[2] /= n;
		}
		double operator[](int i) const { return data[i]; };
		double& operator[](int i) { return data[i]; };
		double data[3];
};
static inline double sqr(double x) {
	return x*x;
}

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator-(const Vector& b) {
	return Vector( -b[0], -b[1], -b[2]);
}
Vector operator*(const double a, const Vector& b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
    
	return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(const Vector &a, const Vector &b)
{
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}
Vector operator/(const Vector& a, const double b) {
	return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

/*---------------------------------------------------------------------------*/
/*                                 Ray class                                 */
/*---------------------------------------------------------------------------*/

class Ray {
	public:
		Vector O, direction;
		Ray(const Vector& O, const Vector& direction) : O(O), direction(direction) {};
};

/*---------------------------------------------------------------------------*/
/*                              Geometry class                               */
/*---------------------------------------------------------------------------*/

class Geometry
{
public:
	Geometry(){};
	virtual bool intersection(const Ray& r, Vector &P, Vector &N, double &t,
                              const Vector& A, const Vector& B, const Vector& C, 
                              double &alpha, double &beta, double &gamma, int &id_mesh) = 0;
	Vector albedo;
	bool isMirror,isTrans;
};

/*---------------------------------------------------------------------------*/
/*                               Sphere class                                */
/*---------------------------------------------------------------------------*/

class Sphere: public Geometry {
	public:
		Vector center;
		double R;

		Sphere(const Vector& center, double R, const Vector& color,bool ism = false, bool ist = false):Geometry(){
			this->center=center;
			this->R=R;
			this->albedo=color;
			this->isMirror=ism;
			this->isTrans=ist;
		}

		bool intersection(const Ray& r, Vector &P, Vector &N, double &t,
								const Vector& A, const Vector& B, const Vector& C, 
								double &alpha, double &beta, double &gamma, int &id_mesh){
			double delta = square(dot(r.direction, r.O-C)) - ((r.O-C).norm2() - R*R);
			if (delta < 0) return false;

			double t1 = dot(r.direction, C - r.O) - sqrt(delta);
			double t2 = dot(r.direction, C - r.O) + sqrt(delta);

			if (t2 < 0) {return false;}
			
			if (t1 > 0) {
				t = t1;
			} else {
				t = t2;
			}

			P = r.O + t * r.direction;
			N = P - center;
			N.normalize();
			return true;
		}
};

/*---------------------------------------------------------------------------*/
/*                               Random Cosine                               */
/*---------------------------------------------------------------------------*/

Vector random_cos(const Vector &N){
    double r1,r2;
    r1 = uniform(engine);
    r2 = uniform(engine);

    double x,y,z;
    x = cos(2*M_PI* r1) * sqrt(1 - r2);
    y = sin(2*M_PI* r1) * sqrt(1 - r2);
    z = sqrt(r2);

    Vector T1,T2;
    double min_num = N[0];
    
	if (N[1]<min_num){min_num = N[1];}
    if (N[2]<min_num){min_num = N[2];}

    if (min_num == N[0]){T1 = Vector(0, N[2], -N[1]);}
    if (min_num == N[1]){T1 = Vector(N[1], 0, -N[0]);}
    if (min_num == N[2]){T1 = Vector(N[1], -N[0], 0);}
    
	T2 = cross(N, T1);
    T1.normalize();
    T2.normalize();

    return x*T1 + y*T2 + z*N;
}

/*---------------------------------------------------------------------------*/
/*                                 Box class                                 */
/*---------------------------------------------------------------------------*/

class BoundingBox {
	public:
		Vector min, max;
		BoundingBox(Vector min = Vector(0, 0, 0), Vector max = Vector(0, 0, 0)) : min(min), max(max) {};
		bool intersect(const Ray &r){
			double t0[3],t1[3];
			double t0_max, t1_min;

			for (int i = 0; i<3; i++){
				double t_min, t_max;
				t_min = std::min( (min[i] - r.O[i]) / r.direction[i], (max[i] - r.O[i]) / r.direction[i] );
				t_max = std::max( (min[i] - r.O[i]) / r.direction[i], (max[i] - r.O[i]) / r.direction[i] );
				t0[i] = t_min;
				t1[i] = t_max;
			}

			t0_max = std::max(t0[0],std::max(t0[1],t0[2]));
			t1_min = std::min(t1[0],std::min(t1[1],t1[2]));

			return (t1_min > t0_max && t0_max > 0);
		}
};

void boxMuller ( double stdev , double &x , double &y ) {
	double r1 = uniform ( engine ) ;
	double r2 = uniform ( engine ) ;
	x = sqrt(-2 * log ( r1 ) ) *cos ( 2 * M_PI*r2 ) *stdev ;
	y = sqrt(-2 * log ( r1 ) ) *sin ( 2 * M_PI*r2 ) *stdev ;
}

/*---------------------------------------------------------------------------*/
/*                              Triangle classes                             */
/*---------------------------------------------------------------------------*/

bool Triangle_Intersection (const Ray& r, const Vector& A, const Vector& B, const Vector& C, double &alpha, double &beta, double &gamma, double &t, Vector &N){
	Vector e1,e2;
	e1 = B-A;
	e2 = C-A;
	N = cross(e1,e2);
	double b1,y1,dv;
	b1 = dot(e2,cross((A - r.O),r.direction));
	y1 = dot(e1,cross((A - r.O),r.direction));
	dv = dot(r.direction,N);

	beta = b1/dv;
	gamma = -y1/dv;
	alpha = 1-beta-gamma;
	t = dot(A-r.O,N)/dot(r.direction,N);

	if (alpha<0 or alpha>1 or beta<0 or beta>1 or gamma<0 or gamma>1 or t<0){
		return false;
	}
	return true;
}

class TriangleIndices {
	public:
		TriangleIndices(int vertex_i = -1, int vertex_j = -1, int vertex_k = -1, int n_i = -1, int n_j = -1, int n_k = -1, int uv_i = -1, int uv_j = -1, int uv_k = -1, int group = -1, bool added = false) : vertex_i(vertex_i), vertex_j(vertex_j), vertex_k(vertex_k), uv_i(uv_i), uv_j(uv_j), uv_k(uv_k), n_i(n_i), n_j(n_j), n_k(n_k), group(group) {};
		int group; 
		int vertex_i, vertex_j, vertex_k;
		int uv_i, uv_j, uv_k;  
		int n_i, n_j, n_k;
};

class BVH{
	public:
		int head,tail;
		BoundingBox box;
		BVH *left, *right;
};

class TriangleMesh: public Geometry {
public:
    ~TriangleMesh() {}
	TriangleMesh(const Vector color) {
        this->albedo = color;
		this->isMirror = false;
		this->isTrans = false;
    };

	bool intersection(const Ray& r, Vector &P, Vector &N, double &t,
                  const Vector& A, const Vector& B, const Vector& C, 
                  double &alpha, double &beta, double &gamma, int &id_mesh) {
		if (!root.box.intersect(r)) {return false;}

		std::list<BVH*> nodes_to_visit;
		nodes_to_visit.push_front(&root);
		while (!nodes_to_visit.empty()) {
			BVH* node = nodes_to_visit.front();
			nodes_to_visit.pop_front();
			if (!node->box.intersect(r)) continue;

			if (node->left == nullptr && node->right == nullptr) {
				for (int i = node->head; i < node->tail; i++) {
					const TriangleIndices& tri = indices[i];
					if (Triangle_Intersection(r, vertices[tri.vertex_i], vertices[tri.vertex_j], vertices[tri.vertex_k], alpha, beta, gamma, t, N)) {
						P = r.O + t * r.direction;
						id_mesh = i;
						std::cout << "Intersection detected with triangle: " << i << std::endl;
						return true;
					}
				}
			} else {
				if (node->left) nodes_to_visit.push_front(node->left);
				if (node->right) nodes_to_visit.push_front(node->right);
			}
		}
		return false;
	}

    void compute_min_max(Vector &min, Vector &max, int l, int r) {
		auto update_min_max = [&](const Vector &v) {
			for (int j = 0; j < 3; ++j) {
				min[j] = std::min(min[j], v[j]);
				max[j] = std::max(max[j], v[j]);
			}
		};

		min = vertices[indices[l].vertex_i];
		max = vertices[indices[l].vertex_i];

		for (int i = l; i < r; ++i) {
			update_min_max(vertices[indices[i].vertex_i]);
			update_min_max(vertices[indices[i].vertex_j]);
			update_min_max(vertices[indices[i].vertex_k]);
    	}
	}


    void recursive_call(BVH *H, int l, int r){
        Vector min,max;
        compute_min_max(min,max,l,r);
        H->box = BoundingBox(min,max);
        H->head = l;
        H->tail = r;
        H->left = NULL;
        H->right = NULL;

        Vector diag = max - min;
        int diag_max = 0;
        if (diag[1] > diag[0]){diag_max = 1;}
        if (diag[2] > diag[diag_max]){diag_max = 2;}

        int pivot_index = l;

        Vector middle = diag*0.5 + min;
        double middle_axis = middle[diag_max];
        for (int i = l; i<r; i++){
            Vector barycenter = (vertices[indices[i].vertex_i] + vertices[indices[i].vertex_j] + vertices[indices[i].vertex_k]) / 3;
            if (barycenter[diag_max] < middle_axis){
                std::swap ( indices [ i ] , indices [ pivot_index ] ) ;
                pivot_index++;
            }
        }

        if (pivot_index-l<=5 or r-pivot_index<=5 or r-l<=10){return;}
        H->left = new BVH;
        H->right = new BVH;
        recursive_call ( H->left , l , pivot_index ) ;
        recursive_call ( H->right , pivot_index, r ) ;
   }

    void factor_move(double s, const Vector& t) {
        for (int i = 0; i < vertices.size(); i++) {
            vertices[i] = vertices[i] * s + t;
        }
    }
	
	void readOBJ(const char* obj) {

		char matfile[255];
		char grp[255];

		FILE* f;
		f = fopen(obj, "r");
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				curGroup++;
			}

			if (line[0] == 'v' && line[1] == ' ') {
				Vector vec;

				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
					col[0] = std::min(1., std::max(0., col[0]));
					col[1] = std::min(1., std::max(0., col[1]));
					col[2] = std::min(1., std::max(0., col[2]));

					vertices.push_back(vec);
					vertexcolors.push_back(col);

				} else {
					sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector vec;
				sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;
				t.group = curGroup;

				char* consumedline = line + 1;
				int offset;

				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vertex_i = vertices.size() + i0; else	t.vertex_i = i0 - 1;
					if (i1 < 0) t.vertex_j = vertices.size() + i1; else	t.vertex_j = i1 - 1;
					if (i2 < 0) t.vertex_k = vertices.size() + i2; else	t.vertex_k = i2 - 1;
					if (j0 < 0) t.uv_i = uvs.size() + j0; else	t.uv_i = j0 - 1;
					if (j1 < 0) t.uv_j = uvs.size() + j1; else	t.uv_j = j1 - 1;
					if (j2 < 0) t.uv_k = uvs.size() + j2; else	t.uv_k = j2 - 1;
					if (k0 < 0) t.n_i = normals.size() + k0; else	t.n_i = k0 - 1;
					if (k1 < 0) t.n_j = normals.size() + k1; else	t.n_j = k1 - 1;
					if (k2 < 0) t.n_k = normals.size() + k2; else	t.n_k = k2 - 1;
					indices.push_back(t);
				} else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vertex_i = vertices.size() + i0; else	t.vertex_i = i0 - 1;
						if (i1 < 0) t.vertex_j = vertices.size() + i1; else	t.vertex_j = i1 - 1;
						if (i2 < 0) t.vertex_k = vertices.size() + i2; else	t.vertex_k = i2 - 1;
						if (j0 < 0) t.uv_i = uvs.size() + j0; else	t.uv_i = j0 - 1;
						if (j1 < 0) t.uv_j = uvs.size() + j1; else	t.uv_j = j1 - 1;
						if (j2 < 0) t.uv_k = uvs.size() + j2; else	t.uv_k = j2 - 1;
						indices.push_back(t);
					} else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vertex_i = vertices.size() + i0; else	t.vertex_i = i0 - 1;
							if (i1 < 0) t.vertex_j = vertices.size() + i1; else	t.vertex_j = i1 - 1;
							if (i2 < 0) t.vertex_k = vertices.size() + i2; else	t.vertex_k = i2 - 1;
							indices.push_back(t);
						} else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vertex_i = vertices.size() + i0; else	t.vertex_i = i0 - 1;
							if (i1 < 0) t.vertex_j = vertices.size() + i1; else	t.vertex_j = i1 - 1;
							if (i2 < 0) t.vertex_k = vertices.size() + i2; else	t.vertex_k = i2 - 1;
							if (k0 < 0) t.n_i = normals.size() + k0; else	t.n_i = k0 - 1;
							if (k1 < 0) t.n_j = normals.size() + k1; else	t.n_j = k1 - 1;
							if (k2 < 0) t.n_k = normals.size() + k2; else	t.n_k = k2 - 1;
							indices.push_back(t);
						}
					}
				}

				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.group = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vertex_i = vertices.size() + i0; else	t2.vertex_i = i0 - 1;
						if (i2 < 0) t2.vertex_j = vertices.size() + i2; else	t2.vertex_j = i2 - 1;
						if (i3 < 0) t2.vertex_k = vertices.size() + i3; else	t2.vertex_k = i3 - 1;
						if (j0 < 0) t2.uv_i = uvs.size() + j0; else	t2.uv_i = j0 - 1;
						if (j2 < 0) t2.uv_j = uvs.size() + j2; else	t2.uv_j = j2 - 1;
						if (j3 < 0) t2.uv_k = uvs.size() + j3; else	t2.uv_k = j3 - 1;
						if (k0 < 0) t2.n_i = normals.size() + k0; else	t2.n_i = k0 - 1;
						if (k2 < 0) t2.n_j = normals.size() + k2; else	t2.n_j = k2 - 1;
						if (k3 < 0) t2.n_k = normals.size() + k3; else	t2.n_k = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					} else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vertex_i = vertices.size() + i0; else	t2.vertex_i = i0 - 1;
							if (i2 < 0) t2.vertex_j = vertices.size() + i2; else	t2.vertex_j = i2 - 1;
							if (i3 < 0) t2.vertex_k = vertices.size() + i3; else	t2.vertex_k = i3 - 1;
							if (j0 < 0) t2.uv_i = uvs.size() + j0; else	t2.uv_i = j0 - 1;
							if (j2 < 0) t2.uv_j = uvs.size() + j2; else	t2.uv_j = j2 - 1;
							if (j3 < 0) t2.uv_k = uvs.size() + j3; else	t2.uv_k = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						} else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vertex_i = vertices.size() + i0; else	t2.vertex_i = i0 - 1;
								if (i2 < 0) t2.vertex_j = vertices.size() + i2; else	t2.vertex_j = i2 - 1;
								if (i3 < 0) t2.vertex_k = vertices.size() + i3; else	t2.vertex_k = i3 - 1;
								if (k0 < 0) t2.n_i = normals.size() + k0; else	t2.n_i = k0 - 1;
								if (k2 < 0) t2.n_j = normals.size() + k2; else	t2.n_j = k2 - 1;
								if (k3 < 0) t2.n_k = normals.size() + k3; else	t2.n_k = k3 - 1;								
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							} else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vertex_i = vertices.size() + i0; else	t2.vertex_i = i0 - 1;
									if (i2 < 0) t2.vertex_j = vertices.size() + i2; else	t2.vertex_j = i2 - 1;
									if (i3 < 0) t2.vertex_k = vertices.size() + i3; else	t2.vertex_k = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								} else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}

			}

		}
		fclose(f);

	}

	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
    BVH root;
	
};

/*---------------------------------------------------------------------------*/
/*                                Scene class                                */
/*---------------------------------------------------------------------------*/

class Scene {
	public:
		Scene() {};
		
		void addObject(Geometry &s){
			objects.push_back(&s);
		}
		
		std::vector<Geometry*> objects;

		bool intersection(const Ray& r, Vector &P, Vector &N, double &t, int &ID, int &id_mesh){
			t = 3E10;
			bool f_inter = false;

			for (int i =0; i<objects.size(); i++){
				Vector local_P,local_N;
				double local_t;
				Vector A, B, C;
				double a, b, c;
				id_mesh = -1;
				int local_id_mesh;

				bool inter = objects[i]->intersection(r, local_P, local_N, local_t, A, B, C, a, b, c, local_id_mesh);
				if (inter and local_t >0){
					f_inter = true;
					if (local_t < t){
						ID = i;
						t = local_t;
						P = local_P;
						N = local_N;
						id_mesh = local_id_mesh;
					}
				}
			}
			return f_inter;
		}

		Vector getColor(const Ray& r, int bounce){
			Vector color(0,0,0);
			Vector P, N;

			if (bounce < 0) return color;

			int objectID;
			double t;
			int id_mesh;

			bool inter = intersection(r, P, N, t, objectID, id_mesh);

			if (inter) {
				if (objects[objectID]->isMirror) {
					Vector R = r.direction - 2 * dot(r.direction, N) * N;
					Ray reflected(P + 0.001 * N, R);
					return getColor(reflected, bounce - 1);
				}

				if (objects[objectID]->isTrans) {
					double n1 = 1;
					double n2 = 1.4;
					Vector correctN = N;
					if (dot(r.direction, N) > 0) {
						std::swap(n1, n2);
						correctN = -correctN;
					}

					double radic = 1 - sqr(n1 / n2) * (1 - sqr(dot(r.direction, correctN)));
					Vector Tt = n1 / n2 * (r.direction - dot(r.direction, correctN) * correctN);
					if (radic < 0) {
						Vector R = r.direction - 2 * dot(r.direction, N) * N;
						Ray reflected(P + 0.001 * N, R);
						return getColor(reflected, bounce - 1);
					}

					Vector Tn = -sqrt(radic) * correctN;
					Vector T = Tt + Tn;
					Ray refracted(P + 0.001 * T, T);
					return getColor(refracted, bounce - 1);
				}

				Vector L(-10, 20, 40);
				double I = 2E10;
				Vector lightVec = (L - P);
				double distlight2 = lightVec.norm2();
				lightVec.normalize();

				Vector Plight, Nlight;
				int objectLight;
				double tlight;
				Ray lightRay(P + 0.001 * N, lightVec);
				double Shadow = 1.0;

				if (intersection(lightRay, Plight, Nlight, tlight, objectLight, id_mesh)) {
					if (tlight * tlight < distlight2) {
						Shadow = 0;
					}
				}

				double intensity = I / (distlight2 * 4 * M_PI);
				Vector L0 = Shadow * intensity * objects[objectID]->albedo / M_PI * std::max(0.0, dot(lightVec, N));
				Ray randomRay(P, random_cos(N));
				color = L0 + objects[objectID]->albedo * getColor(randomRay, bounce - 1);
			}
			return color;
		}
};


int main() {
    auto start = std::chrono::steady_clock::now();
	int W = 512;
	int H = 512;

	double fov = 70 * M_PI / 180;
	Vector Q(0,0,70);

	Sphere Floor   (Vector(0, -1000, 0), 990, Vector(1.0, 0.5, 0.3));
	Sphere Ceiling (Vector(0, 1000, 0),  940, Vector(0.1, 0.8, 0.9));
	Sphere Left    (Vector(-1000, 0, 0), 940, Vector(0.1, 0.8, 0.3));
	Sphere Right   (Vector(1000, 0, 0),  940, Vector(0.5, 0.3, 0.8));
	Sphere Back    (Vector(0, 0, 1000),  900, Vector(0.2, 0.9 ,0.5));
	Sphere Front   (Vector(0, 0, -1000), 940, Vector(0.1, 0.2, 0.6));

	Scene s;
	s.addObject(Floor);
	s.addObject(Ceiling);
	s.addObject(Left);
	s.addObject(Right);
	s.addObject(Back);
	s.addObject(Front);

    TriangleMesh mesh(Vector(0.3,0.3,0.3));
	mesh.readOBJ("cat.obj");
	std::cout << "Number of vertices: " << mesh.vertices.size() << std::endl;
	std::cout << "Number of triangles: " << mesh.indices.size() << std::endl;

	mesh.factor_move(0.6, Vector(0, -10, 0));	
    mesh.recursive_call(&mesh.root, 0 ,mesh.indices.size());

	std::cout << "Transformed vertices:" << std::endl;
	for (const auto& vertex : mesh.vertices) {
		std::cout << vertex[0] << " " << vertex[1] << " " << vertex[2] << std::endl;
	}

    std::cout << mesh.indices.size() << std::endl;
    std::cout << "mesh vertices" << mesh.vertices.size() << std::endl;

    s.addObject(mesh);

	std::vector<unsigned char> image(W * H * 3, 0);
    int nb_paths = 20;
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
            Vector color(0.,0.,0.);
            for (int k = 0; k < nb_paths; k++){
                double x, y;
                boxMuller(1, x, y);
                x = i;
                y = j;
                Vector u(y - W/2 + 0.5, H/2 - x - 0.5, -W / (2 * tan( fov / 2)));
                u.normalize();
                Ray r(Q,u);
                color = color + s.getColor(r,5);   
            }
            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(color[0]/nb_paths, 1. / 2.2));
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(color[1]/nb_paths, 1. / 2.2));
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2]/nb_paths, 1. / 2.2));
        
        }
	}
	stbi_write_png("image8.png", W, H, 3, &image[0], 0);

    auto fin_ish = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(fin_ish - start).count();
    std::cout << "Time is " << elapsed << " microseconds" << std::endl;

	return 0;
}