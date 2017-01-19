
public static class Global {
  public static int N = 20;
}

public static class Util {
  public static int IX(int x, int y) {
    return (x)+(Global.N+2)*(y);
  }

  public static void SWAP(float[] x, float[] x0) {
    for (int a = 0; a < x.length; a++) {
      float temp;
      temp = x[a];
      x[a] = x0[a];
      x0[a] = temp;
    }
  }
}

void add_source ( int N, float[] x, float[] s, float dt )
{
  int i, size=(N+2)*(N+2);
  for ( i=0; i<size; i++ ) x[i] += dt*s[i];
}


void diffuse ( int N, int b, float [] x, float [] x0, float diff, float dt )
{
  int i, j, k;
  float a=dt*diff*N*N;
  for ( k=0; k<20; k++ ) {
    for ( i=1; i<=N; i++ ) {
      for ( j=1; j<=N; j++ ) {
        x[Util.IX(i, j)] = (x0[Util.IX(i, j)] + a*(x[Util.IX(i-1, j)]+x[Util.IX(i+1, j)]
          +x[Util.IX(i, j-1)]+x[Util.IX(i, j+1)]))/(1+4*a);
      }
    }
    //set_bnd ( N, b, x );
  }
}

void advect ( int N, int b, float [] d, float [] d0, float [] u, float [] v, float dt )
{
  int i, j, i0, j0, i1, j1;
  float x, y, s0, t0, s1, t1, dt0;
  dt0 = dt*N;
  for ( i=1; i<=N; i++ ) {
    for ( j=1; j<=N; j++ ) {
      x = i-dt0*u[Util.IX(i, j)]; 
      y = j-dt0*v[Util.IX(i, j)];
      if (x<0.5) x=0.5; 
      if (x>N+0.5) x=N+0.5; 
      i0=(int)x; 
      i1=i0+1; 
      if (y<0.5) y=0.5; 
      if (y>N+0.5) y=N+0.5; 
      j0=(int)y; 
      j1=j0+1; 
      s1 = x-i0; 
      s0 = 1-s1; 
      t1 = y-j0; 
      t0 = 1-t1;
      d[Util.IX(i, j)] = s0*(t0*d0[Util.IX(i0, j0)]+t1*d0[Util.IX(i0, j1)])+
        s1*(t0*d0[Util.IX(i1, j0)]+t1*d0[Util.IX(i1, j1)]);
    }
  }
  //set_bnd ( N, b, d );
}

void dens_step ( int N, float [] x, float [] x0, float [] u, float [] v, float diff, float dt )
{
  add_source ( N, x, x0, dt );
  Util.SWAP ( x0, x ); 
  diffuse ( N, 0, x, x0, diff, dt ); 
  Util.SWAP ( x0, x ); 
  advect ( N, 0, x, x0, u, v, dt );
}

void vel_step ( int N, float [] u, float [] v, float [] u0, float [] v0, 
  float visc, float dt )
{
  add_source ( N, u, u0, dt ); 
  add_source ( N, v, v0, dt );
  Util.SWAP ( u0, u ); 
  diffuse ( N, 1, u, u0, visc, dt );
  Util.SWAP ( v0, v ); 
  diffuse ( N, 2, v, v0, visc, dt );
  project ( N, u, v, u0, v0 );
  Util.SWAP ( u0, u ); 
  Util.SWAP ( v0, v );
  advect ( N, 1, u, u0, u0, v0, dt ); 
  advect ( N, 2, v, v0, u0, v0, dt ); 
  project ( N, u, v, u0, v0 );
}

void project ( int N, float [] u, float [] v, float [] p, float []div )
{
  int i, j, k;
  float h;
  h = 1.0/N;
  for ( i=1; i<=N; i++ ) {
    for ( j=1; j<=N; j++ ) {
      div[Util.IX(i, j)] = -0.5*h*(u[Util.IX(i+1, j)]-u[Util.IX(i-1, j)]+
        v[Util.IX(i, j+1)]-v[Util.IX(i, j-1)]);
      p[Util.IX(i, j)] = 0;
    }
  }

  set_bnd ( N, 0, div ); 
  set_bnd ( N, 0, p );


  for ( k=0; k<20; k++ ) {
    for ( i=1; i<=N; i++ ) {
      for ( j=1; j<=N; j++ ) {
        p[Util.IX(i, j)] = (div[Util.IX(i, j)]+p[Util.IX(i-1, j)]+p[Util.IX(i+1, j)]+p[Util.IX(i, j-1)]+p[Util.IX(i, j+1)])/4 ;
      }
    }
  }
  set_bnd ( N, 0, p ); 


  for (i = 1; i <= N; i++) {
    for (j = 1; j <= N; j++ ) {
      u[Util.IX(i, j)] -= 0.5*(p[Util.IX(i+1, j)]-p[Util.IX(i-1, j)])/h;
      v[Util.IX(i, j)] -= 0.5*(p[Util.IX(i, j+1)]-p[Util.IX(i, j-1)])/h;
    }
  }
  set_bnd ( N, 1, u ); 
  set_bnd ( N, 2, v );
}

void set_bnd ( int N, int b, float [] x )
{
  int i;

  for ( i=1; i<=N; i++ ) {
    x[Util.IX(0, i)] = b==1 ?  -x[Util.IX(1, i)] : x[Util.IX(1, i)]; 
    x[Util.IX(N+1, i)] = b==1 ? -x[Util.IX(N, i)] : x[Util.IX(N, i)]; 
    x[Util.IX(i, 0 )] = b==2 ? -x[Util.IX(i, 1)] : x[Util.IX(i, 1)];
    x[Util.IX(i, N+1)] = b==2 ? -x[Util.IX(i, N)] : x[Util.IX(i, N)];
  }
  x[Util.IX(0, 0 )] = 0.5*(x[Util.IX(1, 0 )]+x[Util.IX(0, 1)]); 
  x[Util.IX(0, N+1)] = 0.5*(x[Util.IX(1, N+1)]+x[Util.IX(0, N )]); 
  x[Util.IX(N+1, 0 )] = 0.5*(x[Util.IX(N, 0 )]+x[Util.IX(N+1, 1)]); 
  x[Util.IX(N+1, N+1)] = 0.5*(x[Util.IX(N, N+1)]+x[Util.IX(N+1, N)]);
}


void setup() {
  size(512, 512);
  frameRate(60);
  background(0);
}