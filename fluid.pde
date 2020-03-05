final int Nx = 168;
final int Ny = 105;
final int iter = 1;
final int SCALE = 4;

int IX(int x, int y) {
  x = constrain(x, 0, Nx-1);
  y = constrain(y, 0, Ny-1);
  return x + (y * Nx);
}

class Fluid {
  int size;
  float dt; //time step
  float diff; //diffusion amount
  float visc; //thickness of fluid

  float[] s; //previous density
  float[] density;

  float[] Vx;
  float[] Vy;

  float[] Vx0;
  float[] Vy0;

  float[] r;
  float[] g;
  float[] b;

  Fluid(float dt, float diffusion, float viscosity) {
    this.dt = dt;
    this.diff = diffusion;
    this.visc = viscosity;

    this.s = new float[Nx*Ny];
    this.density = new float[Nx*Ny];

    this.Vx = new float[Nx*Ny];
    this.Vy = new float[Nx*Ny];

    this.Vx0 = new float[Nx*Ny];
    this.Vy0 = new float[Nx*Ny];

    this.r = new float[Nx*Ny];
    this.g = new float[Nx*Ny];
    this.b = new float[Nx*Ny];
  }

  void step() {
    float visc     = this.visc;
    float diff     = this.diff;
    float dt       = this.dt;
    float[] Vx      = this.Vx;
    float[] Vy      = this.Vy;
    float[] Vx0     = this.Vx0;
    float[] Vy0     = this.Vy0;
    float[] s       = this.s;
    float[] density = this.density;

    diffuse(1, Vx0, Vx, visc, dt);
    diffuse(2, Vy0, Vy, visc, dt);
    
    //project(Vx0, Vy0, Vx, Vy);

    advect(1, Vx, Vx0, Vx0, Vy0, dt);
    advect(2, Vy, Vy0, Vx0, Vy0, dt);

    //project(Vx, Vy, Vx0, Vy0);
    //println(Vx[10000],Vy[10000],Vx0[2500],Vy0[2500],Vx0[400],Vy0[400]);

    diffuse(0, s, r, diff, dt);
    advect(0, r, s, Vx, Vy, dt);
    
    //diffuse(0, s, g, diff, dt);
    //advect(0, g, s, Vx, Vy, dt);
    
    //diffuse(0, s, b, diff, dt);
    //advect(0, b, s, Vx, Vy, dt);
 
  }

  void addDensity(int x, int y, float amount) {
    int index = IX(x, y);
    this.density[index] += amount;
  }

  void addColor(int x, int y, int pick) {
    int index = IX(x, y);
    float dt = this.density[index]-this.s[index];
    /*if (pick==0) {
      this.r[index] += dt;
    }
    if (pick==1) {
      this.g[index] += dt;
    }
    if (pick==2) {
      this.b[index] += dt;
    }*/
    this.r[index] += dt;
  }

  void addVelocity(int x, int y, float amountX, float amountY) {
    int index = IX(x, y);
    this.Vx[index] += amountX;
    this.Vy[index] += amountY;
  }

  void renderD() {
    colorMode(HSB, 255);
    for (int i=0; i<Nx; i++) {
      for (int j=0; j<Ny; j++) {
        float x = i*SCALE;
        float y = j*SCALE;
        float d = this.density[IX(i, j)];
        noStroke();
        fill(d/2, 255, 255);
        square(x, y, SCALE);
      }
    }
  }

  void renderV() {
    for (int i=0; i<Nx; i++) {
      for (int j=0; j<Ny; j++) {
        float x = i*SCALE;
        float y = j*SCALE;
        float vx = this.Vx[IX(i, j)];
        float vy = this.Vy[IX(i, j)];
        stroke(255);
        if (!(abs(vx)<0.1&&abs(vy)<=0.1)) {
          line(x, y, x+vx*SCALE, y+vy*SCALE);
        }
      }
    }
  }

  void fadeD() {
    for (int i=0; i<this.density.length; i++) {
      float d = density[i];
      density[i] = constrain(d-0.1, 0, 255);
    }
  }
}

void diffuse(int b, float[] x, float[] x0, float diff, float dt) {
  float a = dt * diff * (Nx - 2) * (Ny - 2);
  lin_solve(b, x, x0, a, 1 + 6 * a);
}

void lin_solve(int b, float[] x, float[] x0, float a, float c) {
  float cRecip = 1.0 / c;
  for (int k = 0; k < iter; k++) {
    for (int j = 1; j < Ny - 1; j++) {
      for (int i = 1; i < Nx - 1; i++) {
        x[IX(i, j)] =
          (x0[IX(i, j)]
          + a*(    x[IX(i+1, j)]
          +x[IX(i-1, j)]
          +x[IX(i, j+1)]
          +x[IX(i, j-1)]
          )) * cRecip;
      }
    }
    //set_bnd(b, x);
  }
}

void project(float[] velocX, float[] velocY, float[] p, float[] div) {
  for (int j = 1; j < Ny - 1; j++) {
    for (int i = 1; i < Nx - 1; i++) {
      div[IX(i, j)] = -0.5f*(
        (velocX[IX(i+1, j)]-velocX[IX(i-1, j)])/Nx + 
        (velocY[IX(i, j+1)]-velocY[IX(i, j-1)])/Ny
        );
      p[IX(i, j)] = 0;
    }
  }
  set_bnd(0, div); 
  set_bnd(0, p);
  lin_solve(0, p, div, 1, 6);

  for (int j = 1; j < Ny - 1; j++) {
    for (int i = 1; i < Nx - 1; i++) {
      velocX[IX(i, j)] -= 0.5f * (  p[IX(i+1, j)]
        -p[IX(i-1, j)]) * Nx;
      velocY[IX(i, j)] -= 0.5f * (  p[IX(i, j+1)]
        -p[IX(i, j-1)]) * Ny;
    }
  }
  set_bnd(1, velocX);
  set_bnd(2, velocY);
}


void advect(int b, float[] d, float[] d0, float[] velocX, float[] velocY, float dt)
{
  float i0, i1, j0, j1;

  float dtx = dt * (Nx - 2);
  float dty = dt * (Ny - 2);

  float s0, s1, t0, t1;
  float tmp1, tmp2, x, y;

  float Nxfloat = Nx;
  float Nyfloat = Ny;
  float ifloat, jfloat;
  int i, j;

  for (j = 1, jfloat = 1; j < Ny - 1; j++, jfloat++) { 
    for (i = 1, ifloat = 1; i < Nx - 1; i++, ifloat++) {
      tmp1 = dtx * velocX[IX(i, j)];
      tmp2 = dty * velocY[IX(i, j)];
      x    = ifloat - tmp1; 
      y    = jfloat - tmp2;

      if (x < 0.5f) x = 0.5f; 
      if (x > Nxfloat + 0.5f) x = Nxfloat + 0.5f; 
      i0 = floor(x); 
      i1 = i0 + 1.0f;
      if (y < 0.5f) y = 0.5f; 
      if (y > Nyfloat + 0.5f) y = Nyfloat + 0.5f; 
      j0 = floor(y);
      j1 = j0 + 1.0f; 

      s1 = x - i0; 
      s0 = 1.0f - s1; 
      t1 = y - j0; 
      t0 = 1.0f - t1;

      int i0i = int(i0);
      int i1i = int(i1);
      int j0i = int(j0);
      int j1i = int(j1);

      d[IX(i, j)] = s0 * ( t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)])
        + s1 * ( t0 * d0[IX(i1i, j0i)] + t1 * d0[IX(i1i, j1i)]);
    }
  }

  //set_bnd(b, d);
}


void set_bnd(int b, float[] x) {

  for (int i = 1; i < Nx - 1; i++) {
    x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
    x[IX(i, Ny-1)] = b == 2 ? -x[IX(i, Ny-2)] : x[IX(i, Ny-2)];
  }

  for (int j = 1; j < Ny - 1; j++) {
    x[IX(0, j)] = b == 1 ? -x[IX(1, j)] : x[IX(1, j)];
    x[IX(Nx-1, j)] = b == 1 ? -x[IX(Nx-2, j)] : x[IX(Nx-2, j)];
  }

  x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
  x[IX(0, Ny-1)] = 0.5f * (x[IX(1, Ny-1)] + x[IX(0, Ny-2)]);
  x[IX(Nx-1, 0)] = 0.5f * (x[IX(Nx-2, 0)] + x[IX(Nx-1, 1)]);
  x[IX(Nx-1, Ny-1)] = 0.5f * (x[IX(Nx-2, Ny-1)] + x[IX(Nx-1, Ny-2)]);
}
