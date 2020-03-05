
Fluid fluid;

float t = 0;
int pick;

void settings() {
  size(Nx*SCALE, Ny*SCALE);
}

void setup() {
  // Fluid(dt,density,viscosity)
  //dt make fluid get added faster?? not sure
  //viscosity is thickness of fluid
  //if density is used, you do not need fade function
  fluid = new Fluid(0.5, 0.0, 0.00001);
}

void mouseDragged(){
  fluid.addDensity((mouseX/SCALE),(mouseY/SCALE),200);
  fluid.addColor((mouseX/SCALE),(mouseY/SCALE),pick);
  float amtX = mouseX - pmouseX;
  float amtY = mouseY - pmouseY;
  fluid.addVelocity((mouseX/SCALE),(mouseY/SCALE),amtX,amtY);
}

void draw() {
  background(0);
  /*
  int cx = int(0.5*width/SCALE);
  int cy = int(0.5*height/SCALE);
  for (int i=-1; i<=1; i++) {
    for (int j=-1; j<=1; j++) {
      fluid.addDensity(cx+i, cy+j, random(100, 300));
    }
  }
  float angle = noise(t)*TWO_PI*2;
  PVector v = PVector.fromAngle(angle);
  v.mult(1);
  t += 0.01;
  fluid.addVelocity(cx, cy, v.x, v.y);
  */
  
  fluid.step();
  //fluid.renderD();
  //fluid.fadeD();
  /*
  colorMode(HSB, 255);
  for (int i=0; i<Nx; i++) {
    for (int j=0; j<Ny; j++) {
      float x = i*SCALE;
      float y = j*SCALE;
      int index = i + Nx * j;
      float d = fluid.density[index];
      float dt = fluid.density[index]-fluid.s[index];
      noStroke();
      fill(255*dt, d, d);
      square(x, y, SCALE);
    }
  }
  */
  //colorMode(HSB, 255);
  for (int i=0; i<Nx; i++) {
    for (int j=0; j<Ny; j++) {
      float x = i*SCALE;
      float y = j*SCALE;
      int index = i + Nx * j;
      float d = fluid.density[index];
      float r = fluid.r[index];
      float g = fluid.g[index];
      float b = fluid.b[index];
      noStroke();
      fill(r, g, b);
      square(x, y, SCALE);
    }
  }
  
}

void mousePressed(){
  pick = floor(random(3));
}
