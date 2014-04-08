#ifndef I_GAME_H
#define I_GAME_H

namespace Tmpl8 {

#define MAXP1		 256//(65536 >> 2)				// increase to test your optimized code
#define MAXP2		 (MAXP1 << 2)	// because the player is smarter than the AI
#define MAXBULLET	1024

typedef __m128i quint;
typedef __m128 qufl;

class Smoke
{
public:
	struct Puff { 
    union { int x[8]; quint x4[2]; };
    union { int y[8]; quint y4[2]; };
    union { int vy[8]; quint vy4[2]; };
    union { int life[8]; quint life4[2]; };
  };
	Smoke() : frame( 0 ) {};
	void Tick();
  Puff puffs;
	int frame, xpos, ypos; 
};

class Tank
{
public:
	enum { ACTIVE = 1, P1 = 2, P2 = 4 };
	Tank() : pos( float2( 0, 0 ) ), speed( float2( 0, 0 ) ), reloading( 0 ) {};
  ~Tank(){};
	void Fire( unsigned int party, float2& pos, float2& dir );
	void Tick(unsigned int id);
  static float2 targetP1, targetP2;
	float2 pos, speed;
	float maxspeed;
  int flags, reloading, arrayIndex;
};

class Bullet
{
public:
	enum { ACTIVE = 1, P1 = 2, P2 = 4 };
	Bullet() {};
	void Tick(int id);
	float2 pos, speed;
	//int flags;
};

class Surface;
class Surface8;
class Sprite;
class AlignedSprite;
class Game
{
public:
	void SetTarget( Surface* a_Surface ) { m_Surface = a_Surface; }
	void MouseMove( int x, int y ) { m_MouseX = x; m_MouseY = y; }
	void MouseButton( bool b ) { m_LButton = b; }
	void Init();
	void UpdateTanks();
	void UpdateBullets();
	void DrawTanks();
	void PlayerInput();
	void Tick( float a_DT );
  __declspec(align(16)) Surface* m_Surface, *m_Backdrop, *m_Heights, *m_Grid;
	Sprite* m_PXSprite, *m_Smoke;
	int m_ActiveP1, m_ActiveP2;
	int m_MouseX, m_MouseY, m_DStartX, m_DStartY, m_DFrames;
	bool m_LButton, m_PrevButton;
  __declspec(align(16)) Tank* m_Tank;
};

}; // namespace Templ8

#endif