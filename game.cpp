#include "template.h"
#include "game.h"
#include "surface.h"

#define TEST_COLLISION
#define TEST_SMOKE
#define TEST_AIM_AND_SHOOT
#define TEST_MOUNTAINS

using namespace Tmpl8;

// global data (source scope)
static Game* game;

// mountain peaks (push player away)
static float peakx[16] = { 248, 537, 695, 867, 887, 213, 376, 480, 683, 984, 364, 77, 85, 522, 414, 856 };
static float peaky[16] = { 199, 223, 83, 374, 694, 639, 469, 368, 545, 145, 63, 41, 392, 285, 447, 352 };
static float peakh[16] = { 200, 150, 160, 255, 200, 255, 200, 300, 120, 100, 80, 80, 80, 160, 160, 160 };

// player, bullet and smoke data
float2 Tank::targetP1, Tank::targetP2;

float2 pushForces[MAXP1 + MAXP2];

Tank* tankGrid[SCRHEIGHT / 32 + 2][SCRWIDTH / 32 + 2][128];
unsigned int idTankGrid[SCRHEIGHT / 32 + 2][SCRWIDTH / 32 + 2];
int2 tankGridPos[MAXP1 + MAXP2];
int2 tankIpos[MAXP1 + MAXP2];

int tileFlags[SCRHEIGHT / 32 + 2][SCRWIDTH / 32 + 2];

Smoke* smoke[MAXP1 + MAXP2];

static int aliveP1 = MAXP1, aliveP2 = MAXP2;
static Bullet bullet[MAXBULLET];

#ifdef TEST_MOUNTAINS
unsigned long long mountainTiming;
TimerRDTSC mountainTimer;
unsigned long long mountainCount = 0;
#endif
#ifdef TEST_AIM_AND_SHOOT
unsigned long long aimTiming;
TimerRDTSC aimTimer;
unsigned long long aimCount = 0;
#endif
#ifdef TEST_SMOKE
unsigned long long smokeTiming;
TimerRDTSC smokeTimer;
unsigned long long smokeCount = 0;
#endif
// smoke particle effect tick function
void Smoke::Tick()
{
#ifdef TEST_SMOKE
  smokeTimer.Start();
#endif
  //frame is member of Smoke
  unsigned int p = frame >> 3; // Frame divided by 8
  if (frame < 64)
  {
    //if frame % 8 == 0, use p as index. if frame is 56, statement is true and p==7
    if (!(frame++ & 7)) // keep position 
      puff[p].x = xpos,
      puff[p].y = ypos << 8,
      puff[p].vy = -450,
      puff[p].life = 63;
  }

  // First draws 8 puffs, then continues to loop 4 puffs (if i % 2 == 1), so p=1, p=3, p=5, and p=7
  for (unsigned int i = 0; i < p; i++)
  {
    if ((frame < 64) || (i & 1)) // Short-circuiting possible here
    {
      puff[i].x++, puff[i].y += puff[i].vy, puff[i].vy += 3; //Integration of smoke puff

      int animframe = (puff[i].life > 13) ? (9 - (puff[i].life - 14) / 5) : (puff[i].life / 2);
      game->m_Smoke->SetFrame(animframe);
      game->m_Smoke->Draw(puff[i].x - 12, (puff[i].y >> 8) - 12, game->m_Surface);
      if (!--puff[i].life)  // Decrease life, if hits zero, reset smoke and start again
        puff[i].x = xpos,
        puff[i].y = ypos << 8,
        puff[i].vy = -450,
        puff[i].life = 63;
    }
  }

#ifdef TEST_SMOKE
  smokeTimer.Stop();
  smokeTiming += smokeTimer.Interval();
  ++smokeCount;
#endif
}

// bullet Tick function
void Bullet::Tick()
{
  // flags == 0? return
  if (!(flags & Bullet::ACTIVE))
    return;

  // Integration
  float2 prevpos = pos;
  pos += 1.5f * speed, prevpos -= pos - prevpos;

  // Drawing bullet
  game->m_Surface->AddLine(prevpos.x, prevpos.y, pos.x, pos.y, 0x555555);

  // Screen culling
  if ((pos.x < 0) || (pos.x >(SCRWIDTH - 1)) || (pos.y < 0) || (pos.y >(SCRHEIGHT - 1)))
  {
    flags = 0; // off-screen
    return;
  }
  // Determine opponents to check
  const int2 currentTile((int)pos.x >> 5, (int)pos.y >> 5);

  if (!(tileFlags[1 + currentTile.y][1 + currentTile.x] & ((flags & Tank::P1) ? Tank::P2 : Tank::P1)))
    return;

  for (unsigned int i = 0; i < idTankGrid[1+currentTile.y][1+currentTile.x]; i++)
  {
    Tank* t = tankGrid[1 + currentTile.y][1 + currentTile.x][i];
    const int2& tIpos = tankIpos[t->arrayIndex];
    const int2& ipos = *((int2*)&(pos));

    if (t->flags & (Tank::ACTIVE) && t->flags & ((flags & Bullet::P1) ? Tank::P2 : Tank::P1))
    {
      if ((ipos.x >(tIpos.x - 2)) && (ipos.y > (tIpos.y - 2)) &&
        (ipos.x < (tIpos.x + 2)) && (ipos.y < (tIpos.y + 2)))
        continue;

      if (t->flags & Tank::P1) aliveP1--; else aliveP2--; // update counters
      t->flags &= Tank::P1 | Tank::P2;	// kill tank
      smoke[t->arrayIndex]->xpos = tankIpos[t->arrayIndex].x, smoke[t->arrayIndex]->ypos = tankIpos[t->arrayIndex].y;
      flags = 0;						// destroy bullet
      break;
    }
  }
}

// Tank::Fire - spawns a bullet
void Tank::Fire(unsigned int party, float2& pos, float2& dir)
{
  // Low-level: Send tank's flags, saves an OR
  // dir == speed, pos == pos in tank. args unnecessary
  static int currentBullet = 0;
  bullet[currentBullet].flags |= Bullet::ACTIVE + party; // set owner, set active
  bullet[currentBullet].pos = pos, bullet[currentBullet].speed = speed;

  if (++currentBullet == MAXBULLET) currentBullet = 0;
}


// Tank::Tick - update single tank
TimerRDTSC mountainCyclesTimer;
void Tank::Tick(unsigned int id)
{
  if (!(flags & ACTIVE)) // dead tank
  {
    return smoke[id]->Tick();
  }
  
  float2 force = Normalize(((flags & Tank::P1) ? targetP1 : targetP2) - pos);
  
  /*if (tankGridPos[id].x < -1 || tankGridPos[id].x > 33 ||
    tankGridPos[id].y < -1 || tankGridPos[id].y > 25) // in-screen check.
  {
    speed = force;
    pos += speed * (maxspeed * 0.5f);
    return;
  }*/
  
#ifdef TEST_MOUNTAINS
  mountainTimer.Start();
#endif
  // Complexity: O(n*k) (tanks*mountains)
  // evade mountain peaks
  for (unsigned int i = 0; i < 16; i++)
  {
    float2 d(pos.x - peakx[i], pos.y - peaky[i]);
    float sd = (d.x * d.x + d.y * d.y) * 0.2f;
    if (sd < 1500) force += d * 0.03f * (peakh[i] / sd);
  }
#ifdef TEST_MOUNTAINS
  mountainTimer.Stop();
  mountainTiming += mountainTimer.Interval();
  ++mountainCount;
#endif
  //printf("%llu\n", mountainCyclesTimer.Interval());
  // evade other tanks

  force += pushForces[id];

  // High-Level: O(n), it's ok
  // Low-Level: Not really expensive, but has a normalize.
  // SIMD: Nope.
  // evade user dragged line
  if ((flags & P1) && (game->m_LButton))
  {
    float x1 = (float)game->m_DStartX, y1 = (float)game->m_DStartY;
    float x2 = (float)game->m_MouseX, y2 = (float)game->m_MouseY;
    float2 N = Normalize(float2(y2 - y1, x1 - x2));
    float dist = Dot(N, pos) - Dot(N, float2(x1, y1));
    if (fabs(dist) < 10) if (dist > 0) force += 20 * N; else force -= 20 * N;
  }
  // update speed using accumulated force
  float2 dir;
  speed += force, speed = dir = Normalize(speed), pos += speed * (maxspeed * 0.5f);

  tankIpos[id] = int2((int)pos.x, (int)pos.y);

  // shoot, if reloading completed
  if (--reloading >= 0) return;

#ifdef TEST_AIM_AND_SHOOT
  aimTimer.Start();
#endif
  // Calculate possible endpoint on grid
  dir *= 100;
  dir += pos;
  int2 iendpos((int)dir.x >> 5, (int)dir.y >> 5);
  int2 ibegpos = tankGridPos[id];
  
  // Clamp points
  ibegpos.x = MIN(MAX(ibegpos.x, 0), 31);
  ibegpos.y = MIN(MAX(ibegpos.y, 0), 23);
  iendpos.x = MIN(MAX(iendpos.x, 0), 31);
  iendpos.y = MIN(MAX(iendpos.y, 0), 23);
  // Calculate offset
  const int begGtEndX = ibegpos.x > iendpos.x ? -1 : 1;
  const int begGtEndY = ibegpos.y > iendpos.y ? -1 : 1;

  // Check all tiles in between
  for (unsigned int y = ibegpos.y; y != iendpos.y; y += begGtEndY)
  {
    for (unsigned int x = ibegpos.x; x != iendpos.x; x += begGtEndX)
    {
      // If no enemy is present, leave.
      if (!(tileFlags[1 + y][1 + x] & ((flags & Tank::P1) ? Tank::P2 : Tank::P1)))
        continue;

      for (unsigned int i = 0; i < idTankGrid[1 + y][1 + x]; i++)
      {
        if (tankGrid[1 + y][1 + x][i]->flags & ACTIVE
          && tankGrid[1 + y][1 + x][i]->flags & ((flags & Tank::P1) ? Tank::P2 : Tank::P1))
        {
          float2 d = tankGrid[1 + y][1 + x][i]->pos - pos;
          if ((Length(d) < 100) && (Dot(Normalize(d), speed) > 0.99999f))
          {
            Fire(flags & (P1 | P2), pos, speed); // shoot TODO TODO
            reloading = 200; // and wait before next shot is ready
            return;
          }
        }
      }
    }
  }
#ifdef TEST_AIM_AND_SHOOT
  aimTimer.Stop();
  aimTiming += aimTimer.Interval();
  ++aimCount;
#endif
}

int2 glowArrayBounds[40];
// Game::Init - Load data, setup playfield
void Game::Init()
{
  printf("%i\n", sizeof(Tank));
  m_Heights = new Surface("testdata/heightmap.png"), m_Backdrop = new Surface("testdata/backdrop.png"), m_Grid = new Surface(1024, 768);
  Pixel* a1 = m_Grid->GetBuffer(), *a2 = m_Backdrop->GetBuffer(), *a3 = m_Heights->GetBuffer();
  for (int y = 0; y < 768; y++) for (int idx = y * 1024, x = 0; x < 1024; x++, idx++) a1[idx] = (((x & 31) == 0) | ((y & 31) == 0)) ? 0x6600 : 0;
  for (int y = 0; y < 767; y++) for (int idx = y * 1024, x = 0; x < 1023; x++, idx++)
  {
    float3 N = Normalize(float3((float)(a3[idx + 1] & 255) - (a3[idx] & 255), 1.5f, (float)(a3[idx + 1024] & 255) - (a3[idx] & 255))), L(1, 4, 2.5f);
    float h = (float)(a3[x + y * 1024] & 255) * 0.0005f, dx = x - 512.f, dy = y - 384.f, d = sqrtf(dx * dx + dy * dy), dot = Dot(N, Normalize(L));
    int u = max(0, min(1023, (int)(x - dx * h))), v = max(0, min(767, (int)(y - dy * h))), r = (int)Rand(255);
    a2[idx] = AddBlend(a1[u + v * 1024], ScaleColor(ScaleColor(0x33aa11, r) + ScaleColor(0xffff00, (255 - r)), (int)(max(0, dot) * 80.0f) + 10));
  }
  m_Tank = new Tank[MAXP1 + MAXP2];
  m_P1Sprite = new Sprite(new Surface("testdata/p1tank.tga"), 1, Sprite::FLARE);
  m_P2Sprite = new Sprite(new Surface("testdata/p2tank.tga"), 1, Sprite::FLARE);
  m_PXSprite = new Sprite(new Surface("testdata/deadtank.tga"), 1, Sprite::BLACKFLARE);
  m_Smoke = new Sprite(new Surface("testdata/smoke.tga"), 10, Sprite::FLARE);
  // create blue tanks
  Tank::targetP1 = float2(SCRWIDTH / 2, SCRHEIGHT / 2); 
  Tank::targetP2 = float2(SCRWIDTH / 2, SCRHEIGHT / 2); // move to player base
  for (unsigned int i = 0; i < MAXP1; i++)
  {
    Tank& t = m_Tank[i];// = new Tank();
    t.pos = float2((float)((i % 16) * 20), (float)((i / 16) * 20));
    t.speed = float2(0, 0), t.flags = Tank::ACTIVE | Tank::P1, t.maxspeed = (i < (MAXP1 / 2)) ? 0.65f : 0.45f;
    t.arrayIndex = i;
    smoke[t.arrayIndex] = new Smoke();
  }
  // create red tanks
  for (unsigned int i = 0; i < MAXP2; i++)
  {
    Tank& t = m_Tank[i + MAXP1];// = new Tank();
    t.pos = float2((float)((i % 32) * 20 + 700), (float)((i / 32) * 20));
    t.speed = float2(0, 0), t.flags = Tank::ACTIVE | Tank::P2, t.maxspeed = 0.3f;
    t.arrayIndex = MAXP1 + i;
    smoke[t.arrayIndex] = new Smoke();
  }
  game = this; // for global reference
  m_LButton = m_PrevButton = false;
  
  for (int y = 0; y < 40; y++)
  {
    bool firstInRange = false;
    glowArrayBounds[y].y = 40;
    for (int x = 0; x < 40; x++)
    {
      float dist = (x - 20) * (x - 20) + (y - 20) * (y - 20);

      if (dist <= 400) // Sqdst of 20
      {
        if (!firstInRange)
        {
          glowArrayBounds[y].x = x;
          firstInRange = true;
        }
      }
      else if (firstInRange)
      {
        glowArrayBounds[y].y = x;
        break;
      }
    }
  }
}

// Game::DrawTanks - draw the tanks
unsigned long long glowTimings = 0;
unsigned long long restTimings = 0;
unsigned long long glowCount = 0;
void Game::DrawTanks()
{
  static TimerRDTSC timer;

  timer.Start();
  // High-Level: Glow per tank, not per pixel on screen
  // Low-Level: Precalculated & Get Out Early
  for (unsigned int i = 0; i < MAXP1; i++)
  {
    if (m_Tank[i].flags & Tank::ACTIVE)
    {
      //Fixed 4 conversions here by adding tankIpos
      int yMin = tankIpos[i].y - 20;
      int yMax = yMin + 40;
      int yMinOffset = yMin;
      yMin = MAX(yMin, 0);
      yMax = MIN(yMax, SCRHEIGHT);

      int xMin = tankIpos[i].x - 20;
      int xMax = xMin + 40;
      int xMinOffset = xMin;
      xMin = MAX(xMin, 0);
      xMax = MIN(xMax, SCRWIDTH);

      for (int y = yMin; y < yMax; y++)
      {
        const int2& pair = glowArrayBounds[y - yMinOffset];
        int max = xMax;
        if (xMinOffset + pair.y < max)
          max = xMinOffset + pair.y;

        for (int x = xMin + pair.x; x < max; x++)
        {
          m_Surface->GetBuffer()[x + (y << 10)] |= 0x660000;
        }
      }
    }
  }
  timer.Stop();
  glowTimings += timer.Interval();
  ++glowCount;

  timer.Start();
  for (unsigned int i = 0; i < (MAXP1 + MAXP2); i++)
  {
    Tank& t = m_Tank[i];
    float x = t.pos.x, y = t.pos.y;
    const int2& ipos = tankIpos[i];

    if (!(m_Tank[i].flags & Tank::ACTIVE))
      m_PXSprite->Draw(ipos.x - 4, ipos.y - 4, m_Surface); // draw dead tank
    else if (t.flags & Tank::P1) // draw blue tank
    {
      m_P1Sprite->Draw(ipos.x - 4, ipos.y - 4, m_Surface);
      m_Surface->Line(x, y, x + 8 * t.speed.x, y + 8 * t.speed.y, 0x4444ff);
    }
    else // draw red tank
    {
      m_P2Sprite->Draw(ipos.x - 4, ipos.y - 4, m_Surface);
      m_Surface->Line(x, y, x + 8 * t.speed.x, y + 8 * t.speed.y, 0xff4444);
    }

    // Drawing tracks if in screen. 
    if ((x >= 0) && (x < SCRWIDTH) && (y >= 0) && (y < SCRHEIGHT))
      m_Backdrop->GetBuffer()[ipos.x + ipos.y * SCRWIDTH] =
      SubBlend(m_Backdrop->GetBuffer()[ipos.x + ipos.y * SCRWIDTH], 0x030303); // tracks
  }
  timer.Stop();
  restTimings += timer.Interval();

  char glowTimingsStr[128];
  sprintf(glowTimingsStr, "- Glow: %03i", glowTimings / glowCount);
  m_Surface->Print(glowTimingsStr, 20, 70, 0xffff00);

  char drawTimingsStr[128];
  sprintf(drawTimingsStr, "- The rest: %03i", restTimings / glowCount);
  m_Surface->Print(drawTimingsStr, 20, 80, 0xffff00);
}

// Game::PlayerInput - handle player input
void Game::PlayerInput()
{
  if (m_LButton)
  {
    if (!m_PrevButton) m_DStartX = m_MouseX, m_DStartY = m_MouseY, m_DFrames = 0; // start line
    m_Surface->ThickLine(m_DStartX, m_DStartY, m_MouseX, m_MouseY, 0xffffff);
    m_DFrames++;
  }
  else
  {
    if ((m_PrevButton) && (m_DFrames < 150)) // new target location
    Tank::targetP1 = float2((float)m_MouseX, (float)m_MouseY);
    m_Surface->Line(0, (float)m_MouseY, SCRWIDTH - 1, (float)m_MouseY, 0xffffff);
    m_Surface->Line((float)m_MouseX, 0, (float)m_MouseX, SCRHEIGHT - 1, 0xffffff);
  }
  m_PrevButton = m_LButton;
}

//Tank* tankGrid[SCRWIDTH / 32][SCRHEIGHT / 32][128];
//unsigned int idTankGrid[SCRWIDTH / 32][SCRHEIGHT / 32];
unsigned long long utTimingPrev = 0;
unsigned long long blTimingPrev = 0;
unsigned long long dtTimingPrev = 0;
unsigned long long timingCount = 0;

#ifdef TEST_COLLISION
TimerRDTSC collisionTimer;
unsigned long long collisionTime;
unsigned long long colCount; 
#endif
// Game::Tick - main game loop
void Game::Tick(float a_DT)
{
  //printf("%i\n", sizeof(Tank));
  POINT p;
  GetCursorPos(&p);
  ScreenToClient(FindWindow(NULL, "Template"), &p);
  m_LButton = (GetAsyncKeyState(VK_LBUTTON) != 0), m_MouseX = p.x, m_MouseY = p.y;

  TimerRDTSC timer;

  // Drawing backdrop
  memcpy(m_Surface->GetBuffer(), m_Backdrop->GetBuffer(), sizeof(Pixel)* SCRWIDTH * SCRHEIGHT);
  //m_Backdrop->CopyTo(m_Surface, 0, 0);

  // Update tanks
  timer.Start();
  UpdateTanks();
  timer.Stop();
  utTimingPrev += timer.Interval();

  // Update bullets
  timer.Start();
  for (unsigned int i = 0; i < MAXBULLET; i++)
    bullet[i].Tick();
  timer.Stop();
   blTimingPrev += timer.Interval();

  // Draw tanks
  timer.Start();
  DrawTanks();
  timer.Stop();
  dtTimingPrev += timer.Interval();

  PlayerInput();

  //Timings:
  timingCount++;

  char utTimingsStr[128];
  sprintf(utTimingsStr, "Update Tanks: %03i", utTimingPrev / timingCount);
  m_Surface->Print(utTimingsStr, 10, 40, 0xffff00);

  char blTimingsStr[128];
  sprintf(blTimingsStr, "Update Bullets: %03i", blTimingPrev / timingCount);
  m_Surface->Print(blTimingsStr, 10, 50, 0xffff00);

  char dtTimingsStr[128];
  sprintf(dtTimingsStr, "Draw Tanks: %03i", dtTimingPrev / timingCount);
  m_Surface->Print(dtTimingsStr, 10, 60, 0xffff00);

#ifdef TEST_SMOKE
  if (smokeCount > 0)
  {
    char smokeTimingstr[128];
    sprintf(smokeTimingstr, "Smoke::Tick(): %03i", smokeTiming / smokeCount);
    game->m_Surface->Print(smokeTimingstr, 10, 90, 0xffff00);
  }
#endif

#ifdef TEST_COLLISION
  char tankColstr[128];
  sprintf(tankColstr, "Tank-Tank Col.: %0llu", collisionTime / ++colCount);
  game->m_Surface->Print(tankColstr, 10, 100, 0xffff00);
#endif
#ifdef TEST_AIM_AND_SHOOT
  char aimstr[128];
  sprintf(aimstr, "Aim and shoot: %0llu", aimTiming / aimCount);
  game->m_Surface->Print(aimstr, 10, 110, 0xffff00);
#endif
#ifdef TEST_MOUNTAINS
  char mtstr[128];
  sprintf(mtstr, "Mountain push: %0llu", mountainTiming / mountainCount);
  game->m_Surface->Print(mtstr, 10, 120, 0xffff00);
#endif

  char buffer[128];
  if ((aliveP1 > 0) && (aliveP2 > 0))
  {
    sprintf(buffer, "blue army: %03i  red army: %03i", aliveP1, aliveP2);
    return m_Surface->Print(buffer, 10, 10, 0xffff00);
  }
  if (aliveP1 == 0)
  {
    sprintf(buffer, "sad, you lose... red left: %i", aliveP2);
    return m_Surface->Print(buffer, 200, 370, 0xffff00);
  }
  sprintf(buffer, "nice, you win! blue left: %i", aliveP1);
  m_Surface->Print(buffer, 200, 370, 0xffff00);

}

// Stolen from http://en.wikipedia.org/wiki/Fast_inverse_square_root
float Q_rsqrt(float number)
{
  long i;
  float x2, y;
  const float threehalfs = 1.5F;

  x2 = number * 0.5F;
  y = number;
  i = *(long *)&y;                       // evil floating point bit level hacking
  i = 0x5f3759df - (i >> 1);               // what the fuck?
  y = *(float *)&i;
  y = y * (threehalfs - (x2 * y * y));   // 1st iteration
  //      y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed

  return y;
}

void Game::UpdateTanks()
{
  // Clear array
  memset(pushForces, 0, sizeof(pushForces));
  memset(idTankGrid, 0, sizeof(idTankGrid));
  memset(tileFlags, 0, sizeof(tileFlags));
  
  //memset(tankGrid, 0, sizeof(Tank*)* ((SCRWIDTH / 32 + 2) * (SCRHEIGHT / 32 + 2) * 128));
  for (unsigned int i = 0; i < (MAXP1 + MAXP2); i++)
  {
    tankGridPos[i] = int2(m_Tank[i].pos.x, m_Tank[i].pos.y);
    int2& ipos = tankGridPos[i];
    ipos.x >>= 5;
    ipos.y >>= 5;

    if (ipos.x >= -1 && ipos.x < 33 && ipos.y >= -1 && ipos.y < 25)
    {
      tankGrid[1 + ipos.y][1 + ipos.x][idTankGrid[1 + ipos.y][1 + ipos.x]++] = &game->m_Tank[i];
      tileFlags[1 + ipos.y][1 + ipos.x] |= game->m_Tank[i].flags;
    }
  }

#ifdef TEST_COLLISION
  collisionTimer.Start();
#endif

  for (unsigned int i = 0; i < (MAXP1 + MAXP2); i++)
  {
    if (!(m_Tank[i].flags & Tank::ACTIVE))
      continue;

    const int2& ipos = tankGridPos[i];
    if (ipos.x < -1 || ipos.x > 32 || ipos.y < -1 || ipos.y > 24)
      continue;

    // Current grid tile
    for (int j = 0; j < idTankGrid[1 + ipos.y][1 + ipos.x]; j++)
    {
      if (&m_Tank[i] == tankGrid[1 + ipos.y][1 + ipos.x][j])
        continue;

      float2 d = m_Tank[i].pos - tankGrid[1 + ipos.y][1 + ipos.x][j]->pos;
      float len = Dot(d, d);// Length(d);
      if (len < 8 * 8)
      {
        len = Q_rsqrt(len);// Inv. sqrt
        auto dnorm = d * len * 2.0f;
        pushForces[i] += dnorm;
        pushForces[tankGrid[1 + ipos.y][1 + ipos.x][j]->arrayIndex] -= dnorm;
      }
      else if (len < 16 * 16)
      {
        len = Q_rsqrt(len);// Inv. sqrt
        auto dnorm = d * len * 0.4f;
        pushForces[i] += dnorm;
        pushForces[tankGrid[1 + ipos.y][1 + ipos.x][j]->arrayIndex] -= dnorm;
      }
    }
    if (ipos.x > -1)
    {
      // Left grid tile
      for (int j = 0; j < idTankGrid[1 + ipos.y][1 + ipos.x - 1]; j++)
      {
        float2 d = m_Tank[i].pos - tankGrid[1 + ipos.y][1 + ipos.x - 1][j]->pos;
        float len = Dot(d,d);
        if (len < 8 * 8)
        {
          len = Q_rsqrt(len);// Inv. sqrt
          auto dnorm = d * len * 2.0f;
          pushForces[i] += dnorm;
          pushForces[tankGrid[1 + ipos.y][1 + ipos.x - 1][j]->arrayIndex] -= dnorm;
        }
        else if (len < 16 * 16)
        {
          len = Q_rsqrt(len);// Inv. sqrt
          auto dnorm = d * len * 0.4f;
          pushForces[i] += dnorm;
          pushForces[tankGrid[1 + ipos.y][1 + ipos.x - 1][j]->arrayIndex] -= dnorm;
        }
      }
      // Bottom-Left grid tile
      /*if (ipos.y < 24)
      {
        for (int j = 0; j < idTankGrid[1 + ipos.y + 1][1 + ipos.x-1]; j++)
        {
          float2 d = m_Tank[i].pos - tankGrid[1 + ipos.y + 1][1 + ipos.x-1][j]->pos;
          float len = Dot(d,d);
          if (len < 8 * 8)
          {
            len = Q_rsqrt(len);// Inv. sqrt
            auto dnorm = d * len * 2.0f;
            pushForces[i] += dnorm;
            pushForces[tankGrid[1 + ipos.y+1][1 + ipos.x - 1][j]->arrayIndex] -= dnorm;
          }
          else if (len < 16 * 16)
          {
            len = Q_rsqrt(len);// Inv. sqrt
            auto dnorm = d * len * 0.4f;
            pushForces[i] += dnorm;
            pushForces[tankGrid[1 + ipos.y+1][1 + ipos.x - 1][j]->arrayIndex] -= dnorm;
          }
        }
      }*/
    }
    if (ipos.x < 32)
    {
      // Right grid tile
      for (int j = 0; j < idTankGrid[1 + ipos.y][1 + ipos.x + 1]; j++)
      {
        float2 d = m_Tank[i].pos - tankGrid[1 + ipos.y][1 + ipos.x + 1][j]->pos;
        float len = Dot(d,d);
        if (len < 8 * 8)
        {
          len = Q_rsqrt(len);// Inv. sqrt
          auto dnorm = d * len * 2.0f;
          pushForces[i] += dnorm;
          pushForces[tankGrid[1 + ipos.y][1 + ipos.x + 1][j]->arrayIndex] -= dnorm;
        }
        else if (len < 16 * 16)
        {
          len = Q_rsqrt(len);// Inv. sqrt
          auto dnorm = d * len * 0.4f;
          pushForces[i] += dnorm;
          pushForces[tankGrid[1 + ipos.y][1 + ipos.x + 1][j]->arrayIndex] -= dnorm;
        }
      }
      // Bottom-Right grid tile
      if (ipos.y < 24)
      {
        for (int j = 0; j < idTankGrid[1 + ipos.y + 1][1 + ipos.x + 1]; j++)
        {
          float2 d = m_Tank[i].pos - tankGrid[1 + ipos.y + 1][1 + ipos.x + 1][j]->pos;
          float len = Dot(d,d);
          if (len < 8 * 8)
          {
            len = Q_rsqrt(len);// Inv. sqrt
            auto dnorm = d * len * 2.0f;
            pushForces[i] += dnorm;
            pushForces[tankGrid[1 + ipos.y + 1][1 + ipos.x + 1][j]->arrayIndex] -= dnorm;
          }
          else if (len < 16 * 16)
          {
            len = Q_rsqrt(len);// Inv. sqrt
            auto dnorm = d * len * 0.4f;
            pushForces[i] += dnorm;
            pushForces[tankGrid[1 + ipos.y + 1][1 + ipos.x + 1][j]->arrayIndex] -= dnorm;
          }
        }
      }
    }
    // Bottom grid tile
    if (ipos.y < 24)
    {
      for (int j = 0; j < idTankGrid[1 + ipos.y + 1][1 + ipos.x]; j++)
      {
        float2 d = m_Tank[i].pos - tankGrid[1 + ipos.y + 1][1 + ipos.x][j]->pos;
        float len = Dot(d,d);
        if (len < 8 * 8)
        {
          len = Q_rsqrt(len);// Inv. sqrt
          auto dnorm = d * len * 2.0f;
          pushForces[i] += dnorm;
          pushForces[tankGrid[1 + ipos.y + 1][1 + ipos.x][j]->arrayIndex] -= dnorm;
        }
        else if (len < 16 * 16)
        {
          len = Q_rsqrt(len);// Inv. sqrt
          auto dnorm = d * len * 0.4f;
          pushForces[i] += dnorm;
          pushForces[tankGrid[1 + ipos.y + 1][1 + ipos.x][j]->arrayIndex] -= dnorm;
        }
      }
    }
    // Bottom grid tile
    if (ipos.y > -1)
    {
      for (int j = 0; j < idTankGrid[1 + ipos.y - 1][1 + ipos.x]; j++)
      {
        float2 d = m_Tank[i].pos - tankGrid[1 + ipos.y - 1][1 + ipos.x][j]->pos;
        float len = Dot(d,d);
        if (len < 8 * 8)
        {
          len = Q_rsqrt(len);// Inv. sqrt
          auto dnorm = d * len * 2.0f;
          pushForces[i] += dnorm;
          pushForces[tankGrid[1 + ipos.y - 1][1 + ipos.x][j]->arrayIndex] -= dnorm;
        }
        else if (len < 16 * 16)
        {
          len = Q_rsqrt(len);// Inv. sqrt
          auto dnorm = d * len * 0.4f;
          pushForces[i] += dnorm;
          pushForces[tankGrid[1 + ipos.y - 1][1 + ipos.x][j]->arrayIndex] -= dnorm;
        }
      }
    }
  }
#ifdef TEST_COLLISION
  collisionTimer.Stop();
  collisionTime += collisionTimer.Interval();
#endif
  for (unsigned int i = 0; i < (MAXP1 + MAXP2); i++)
    m_Tank[i].Tick(i);
}
