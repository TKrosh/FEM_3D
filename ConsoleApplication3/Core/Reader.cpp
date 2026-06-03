#include "Reader.h"

void Reader::ReadMeshFromFile(const char* FileName, Mesh& TaskMesh)
{
	std::ifstream InfoFile;
	InfoFile.open(FileName);

	int n, material, Bindex;
	// time to start read Points
	InfoFile >> n;
	TaskMesh.SetNumberOfVertexes(n);
	for (int i = 0; i < n; i++)
	{
		double x, y, z;
		InfoFile >> x >> y >> z;
		TaskMesh.AddVertex(x, y, z);
	}

	// time to start read elements
	std::vector<int> indexes{0, 0, 0, 0, 0, 0, 0, 0};
	InfoFile >> n;
	TaskMesh.SetNumberOfElements(n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < 8; j++) InfoFile >> indexes[j];
		InfoFile >> material;

		/* we need some way to find : hx, hy, hz;
		 let's do that little bit simple
		if it doesn't work will see
		 
			7------6
           /      /|
		  /      / |
		 4------5  |
		 |  3---|--2
	   hx| / hy | /        z  y
		 |/ hx  |/         |/
		 0------1          *-> x
		 */
		Vector3D P0 = TaskMesh.Vertexes[indexes[0]];
		double hx = TaskMesh.Vertexes[indexes[1]].X - P0.X;
		double hy = TaskMesh.Vertexes[indexes[3]].Y - P0.Y;
		double hz = TaskMesh.Vertexes[indexes[4]].Z - P0.Z;

		TaskMesh.AddElement(indexes, material, hx, hy, hz);
	}

	// time to start read boundary elements
	std::vector<int> Boundaryindexes{ 0, 0, 0, 0};
	InfoFile >> n;
	TaskMesh.SetNumberOfBoundaryElements(n);
	for (int i = 0; i < n; i++)
	{
		/* we need some way to find : hx, hy;
		 let's do that little bit simple
		 if it doesn't work will see

		  e2
		  3------2
		  |      |
		  |      |
		  |      |
		  0------1 e1
		*/

		for (int j = 0; j < 4; j++) InfoFile >> Boundaryindexes[j];
		InfoFile >> material >> Bindex;

		Vector3D P0 = TaskMesh.Vertexes[Boundaryindexes[0]];
		Vector3D e1 = TaskMesh.Vertexes[Boundaryindexes[1]] - P0;
		Vector3D e2 = TaskMesh.Vertexes[Boundaryindexes[3]] - P0;

		TaskMesh.AddBoundaryElement(Boundaryindexes, material, Bindex, e1, e2);
	}
}

void Reader::BuildcubeMesh(const char* FileName, int n)
{
    std::ofstream out;
    out.open(FileName);
    const int N = n + 1;               // points per side
    const int layerSize = N * N;       // points in one z-layer
    double inv_n = 1.0 / n;   // или static_cast<double>(n)

    std::vector<std::vector<int>> layerIdx(N, std::vector<int>(N, -1));
    std::vector<std::pair<int, int>> localPoints(layerSize);
    int idx = 0;
    for (int r = 0; 2 * r < N; ++r) {
        int minXY = r, maxXY = n - r;
        if (minXY == maxXY) {
            layerIdx[minXY][minXY] = idx;
            localPoints[idx++] = { minXY, minXY };
        }
        else {
            layerIdx[minXY][minXY] = idx;
            localPoints[idx++] = { minXY, minXY };
            for (int x = minXY + 1; x <= maxXY; ++x) {
                layerIdx[x][minXY] = idx;
                localPoints[idx++] = { x, minXY };
            }
            for (int y = minXY + 1; y <= maxXY; ++y) {
                layerIdx[maxXY][y] = idx;
                localPoints[idx++] = { maxXY, y };
            }
            for (int x = maxXY - 1; x >= minXY; --x) {
                layerIdx[x][maxXY] = idx;
                localPoints[idx++] = { x, maxXY };
            }
            for (int y = maxXY - 1; y >= minXY + 1; --y) {
                layerIdx[minXY][y] = idx;
                localPoints[idx++] = { minXY, y };
            }
        }
    }

    // ----- output vertices -----
    out << (N * N * N) << "\n";
    for (int z = 0; z <= n; ++z) {
        for (const auto& p : localPoints) {
            out << (p.first * inv_n) << " " << (p.second * inv_n) << " " << (z * inv_n) << "\n";
        }
    }

    // ----- output cubes (8 vertex indices each) -----
    out << (n * n * n) << "\n";
    for (int k = 0; k < n; ++k) {
        int base0 = k * layerSize;
        int base1 = (k + 1) * layerSize;
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                int v1 = base0 + layerIdx[i][j];
                int v2 = base0 + layerIdx[i + 1][j];
                int v3 = base0 + layerIdx[i + 1][j + 1];
                int v4 = base0 + layerIdx[i][j + 1];
                int v5 = base1 + layerIdx[i][j];
                int v6 = base1 + layerIdx[i + 1][j];
                int v7 = base1 + layerIdx[i + 1][j + 1];
                int v8 = base1 + layerIdx[i][j + 1];
                out << v1 << " " << v2 << " " << v3 << " " << v4 << " "
                    << v5 << " " << v6 << " " << v7 << " " << v8 << " " << 0 << "\n";
            }
        }
    }

    // ----- output external lateral faces (6 * n * n quads) -----
    out << (6 * n * n) << "\n";
    auto outFace = [&](int a, int b, int c, int d) {
        out << a << " " << b << " " << c << " " << d << " " << 0 << " " << 0 << "\n";
        };

    // -y (y = 0)
    for (int z = 0; z < n; ++z) {
        for (int x = 0; x < n; ++x) {
            outFace(z * layerSize + layerIdx[x][0],
                z * layerSize + layerIdx[x + 1][0],
                (z + 1) * layerSize + layerIdx[x + 1][0],
                (z + 1) * layerSize + layerIdx[x][0]);
        }
    }

    // +x (x = n)
    for (int z = 0; z < n; ++z) {
        for (int y = 0; y < n; ++y) {
            outFace(z * layerSize + layerIdx[n][y],
                z * layerSize + layerIdx[n][y + 1],
                (z + 1) * layerSize + layerIdx[n][y + 1],
                (z + 1) * layerSize + layerIdx[n][y]);
        }
    }

    // +y (y = n)
    for (int z = 0; z < n; ++z) {
        for (int x = 0; x < n; ++x) {
            outFace(z * layerSize + layerIdx[x + 1][n],
                z * layerSize + layerIdx[x][n],
                (z + 1) * layerSize + layerIdx[x][n],
                (z + 1) * layerSize + layerIdx[x + 1][n]);
        }
    }

    // -x (x = 0)
    for (int z = 0; z < n; ++z) {
        for (int y = 0; y < n; ++y) {
            outFace(z * layerSize + layerIdx[0][y + 1],
                z * layerSize + layerIdx[0][y],
                (z + 1) * layerSize + layerIdx[0][y],
                (z + 1) * layerSize + layerIdx[0][y + 1]);
        }
    }

    // +z (z = n)
    {
        int baseZ = n * layerSize;
        for (int y = 0; y < n; ++y) {
            for (int x = 0; x < n; ++x) {
                outFace(baseZ + layerIdx[x][y],
                    baseZ + layerIdx[x + 1][y],
                    baseZ + layerIdx[x + 1][y + 1],
                    baseZ + layerIdx[x][y + 1]);
            }
        }
    }

    // -z (z = 0)
    {
        int baseZ = 0;
        for (int y = 0; y < n; ++y) {
            for (int x = 0; x < n; ++x) {
                outFace(baseZ + layerIdx[x][y],
                    baseZ + layerIdx[x][y + 1],
                    baseZ + layerIdx[x + 1][y + 1],
                    baseZ + layerIdx[x + 1][y]);
            }
        }
    }

    out.close();
    std::cout << "Complete Mesh saved in file:" << FileName;
}
