/*
This code computes the Voronoi diagram of a set of points with weights.
In the case of any code similarity, this code was implemented diectly during the live coding session of the course.
I hope this clears up any confusion in terms of similarity in code.
*/

#include <vector>
#include <string>

/*
================================= Vector Class =================================
*/
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

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
};
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
};
Vector operator*(const double a, const Vector& b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
};
Vector operator*(const Vector& a, const double b) {
	return Vector(a[0]*b, a[1]*b, a[2]*b);
};
Vector operator/(const Vector& a, const double b) {
	return Vector(a[0] / b, a[1] / b, a[2] / b);
};
double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
};
Vector cross(const Vector& a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
};
Vector operator-(const Vector& a) {
	return Vector(-a[0], -a[1], -a[2]);
}

class Polygon {  
public:
	std::vector<Vector> vertices;
};	

/*
================================= Saving images =================================
This function is used to save the Voronoi diagram as an SVG file.
NB: This function came from external source.
*/
void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none") {
	FILE* f = fopen(filename.c_str(), "w+"); 
	fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
	for (int i=0; i<polygons.size(); i++) {
		fprintf(f, "<g>\n");
		fprintf(f, "<polygon points = \""); 
		for (int j = 0; j < polygons[i].vertices.size(); j++) {
			fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
		}
		fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
		fprintf(f, "</g>\n");
	}
	fprintf(f, "</svg>\n");
	fclose(f);
};

/*
================================= Clipping function =================================
`clip_by_bisector` uses the Sutherland-Hodgman algorithm to clip a cell by a bisector of UV.
*/
Polygon clip_by_bisector(const Polygon &cell, const Vector& P0, double w0, const Vector& Pi, double wi) {
	Polygon result;
	Vector offset = (w0 - wi)/(2.*(P0 - Pi).norm2()) * (Pi - P0);
	Vector M = (P0 + Pi) * 0.5;
	Vector Mprime = offset + M;
	int N = cell.vertices.size();

	for (int i = 0; i < N; i++) {
		const Vector& A = cell.vertices[i==0 ? (N-1) : i-1];
		const Vector& B = cell.vertices[i];
		if ((B-P0).norm2() - w0 <= (B-Pi).norm2() - wi) { // B is inside the voronoi of P0
			if ((A-P0).norm2() - w0 > (A-Pi).norm2() - wi) { // A is outside
				double t = dot(Mprime - A, Pi - P0) / dot(B - A, Pi - P0);
				Vector P = A + t * (B - A);
				result.vertices.push_back(P);
			}
			result.vertices.push_back(B);
		} else {
			if ((A-P0).norm2() -w0 <= (A-Pi).norm2() - wi) { // A is inside
				double t = dot(Mprime - A, Pi - P0) / dot(B - A, Pi - P0);
				Vector P = A + t * (B - A);
				result.vertices.push_back(P);
			}
		}
	}
	return result;
};

/*
================================= Voronoi Class =================================
This class represents the Voronoi diagram of a set of points in 2D space.
*/
class Voronoi {
public:
	Voronoi(){ };

	void compute() {
		Polygon square;
		square.vertices.push_back(Vector(0, 0, 0));
		square.vertices.push_back(Vector(1, 0, 0));
		square.vertices.push_back(Vector(1, 1, 0));
		square.vertices.push_back(Vector(0, 1, 0));

		cells.resize(points.size());
		for (int i = 0; i < points.size(); i++) {
			Polygon cell = square;
			for (int j = 0; j < points.size(); j++) {
				if (i == j) continue;
				cell = clip_by_bisector(cell, points[i], weights[i], points[j], weights[j]);
			}
			cells[i] = cell;
		}
	};

	std::vector<Vector> points;
	std::vector<Polygon> cells;
	std::vector<double> weights;
};

int main() {
	Voronoi vor;
	int N = 10000;
	vor.points.resize(N);
	vor.weights.resize(N);
	for (int i=0; i < N; i++) {
		vor.points[i] = Vector(rand()/double(RAND_MAX), rand()/double(RAND_MAX), 0);
		vor.weights[i] = rand()/double(RAND_MAX);
	}
	vor.compute();

	save_svg(vor.cells, "vor.svg");
	return 0;
}
