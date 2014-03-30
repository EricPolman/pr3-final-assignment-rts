#include "template.h"
#include "game.h"
#include "surface.h"

using namespace Tmpl8;

// global data (source scope)
static Game* game;

// mountain peaks (push player away)
static float peakx[16] = { 248, 537, 695, 867, 887, 213, 376, 480, 683, 984, 364, 77, 85, 522, 414, 856 };
static float peaky[16] = { 199, 223, 83, 374, 694, 639, 469, 368, 545, 145, 63, 41, 392, 285, 447, 352 };
static float peakh[16] = { 200, 150, 160, 255, 200, 255, 200, 300, 120, 100, 80, 80, 80, 160, 160, 160 };

// player, bullet and smoke data
static int aliveP1 = MAXP1, aliveP2 = MAXP2;
static Bullet bullet[MAXBULLET];

static const int qtMaxDepth = 2;
static const int qtMaxTanks = 16;
struct QuadTree
{
  QuadTree* nodes[4];
  Tank* tanks[MAXP1 + MAXP2];
  const int boundXmin, boundXmax, boundYmin, boundYmax;
  int currentAmountOfTanks;
  int currentDepth;

  QuadTree(int bXmin, int bXmax, int bYmin, int bYmax, int depth) :
    boundXmin(bXmin), boundXmax(bXmax), boundYmin(bYmin), boundYmax(bYmax),
    currentDepth(depth), currentAmountOfTanks(0){
    nodes[0] = nullptr;
  };

  ~QuadTree()
  {
    if (nodes[0] != nullptr)
    {
      delete nodes[0];
      delete nodes[1];
      delete nodes[2];
      delete nodes[3];
    }
  }

  void Split()
  {
    nodes[0] = new QuadTree(boundXmin, boundXmax >> 1, boundYmin, boundYmax >> 1, currentDepth + 1);
    nodes[1] = new QuadTree(boundXmax >> 1, boundXmax, boundYmin, boundYmax >> 1, currentDepth + 1);
    nodes[2] = new QuadTree(boundXmin, boundXmax >> 1, boundYmax >> 1, boundYmax, currentDepth + 1);
    nodes[3] = new QuadTree(boundXmax >> 1, boundXmax, boundYmax >> 1, boundYmax, currentDepth + 1);

    for (int i = 0; i < currentAmountOfTanks; i++)
    {
      nodes[0]->Add(tanks[i]);
      nodes[1]->Add(tanks[i]);
      nodes[2]->Add(tanks[i]);
      nodes[3]->Add(tanks[i]);
    }
  }
  void Add(Tank* a_tank)
  {
    if (Intersects(a_tank))
    {
      if (currentAmountOfTanks < qtMaxTanks)
      {
        tanks[currentAmountOfTanks++] = a_tank;
      }
      else
      {
        if (nodes[0] == nullptr)
        {
          if (currentDepth < qtMaxDepth)
            Split();
          else
            tanks[currentAmountOfTanks++] = a_tank;
        }
        else
        {
          nodes[0]->Add(a_tank);
          nodes[1]->Add(a_tank);
          nodes[2]->Add(a_tank);
          nodes[3]->Add(a_tank);
        }
      }
    }
  }
  bool Intersects(Tank* a_tank)
  {
    float tankX = a_tank->pos.x, tankY = a_tank->pos.y;
    if (tankX - 16 > boundXmax || tankX + 16 < boundXmin ||
      tankY - 16 > boundYmax || tankY + 16 < boundYmin)
      return false;
    return true;
  }

  void Collide(Tank* a_tank, float2& a_force)
  {
    if (Intersects(a_tank))
    {
      if (nodes[0] == nullptr)
      {
        for (int i = 0; i < currentAmountOfTanks; i++)
        {
          if (tanks[i] == a_tank) continue;
          float2 d = a_tank->pos - tanks[i]->pos;
          if (Length(d) < 8) a_force += Normalize(d) * 2.0f;
          else if (Length(d) < 16) a_force += Normalize(d) * 0.4f;
        }
      }
      else
      {
        nodes[0]->Collide(a_tank, a_force);
        nodes[1]->Collide(a_tank, a_force);
        nodes[2]->Collide(a_tank, a_force);
        nodes[3]->Collide(a_tank, a_force);
      }
    }
  }
};
QuadTree* quadTree;

// smoke particle effect tick function
void Smoke::Tick()
{
  //frame is member of Smoke
  unsigned int p = frame >> 3; // Frame divided by 8
  if (frame < 64)
  {
    //if frame % 8 == 0, use p as index. if frame is 56, statement is true and p==7
    if (!(frame++ & 7)) //Clamp frame to range of 8 puffs, check if 0, increment value anyway
      puff[p].x = xpos,
      puff[p].y = ypos << 8,
      puff[p].vy = -450,
      puff[p].life = 63;
  }

  // First draws 8 puffs, then continues to loop 4 puffs (if i % 2 == 1), so p=1, p=3, p=5, and p=7
  for (unsigned int i = 0; i < p; i++)
  {
    if ((frame < 64) || (i & 1))
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
    flags = 0; // off-screen

  // Determine opponents to check
  unsigned int start = 0, end = MAXP1;
  if (flags & P1) start = MAXP1, end = MAXP1 + MAXP2;

  // Distance checking?
  for (unsigned int i = start; i < end; i++) // check all opponents
  {
    Tank* t = game->m_Tank[i];

    if (!((t->flags & Tank::ACTIVE) &&
      (pos.x >(t->pos.x - 2)) && (pos.y > (t->pos.y - 2)) &&
      (pos.x < (t->pos.x + 2)) && (pos.y < (t->pos.y + 2))))
      continue;

    if (t->flags & Tank::P1) aliveP1--; else aliveP2--; // update counters
    t->flags &= Tank::P1 | Tank::P2;	// kill tank
    flags = 0;						// destroy bullet
    break;
  }
}

// Tank::Fire - spawns a bullet
void Tank::Fire(unsigned int party, float2& pos, float2& dir)
{
  // High-level: Keep track of next bullet position. O(n*k) to O(n)
  // Low-level: Send tank's flags, saves an OR
  // dir == speed, pos == pos in tank. args unnecessary
  static int currentBullet = 0;
  bullet[currentBullet].flags |= Bullet::ACTIVE + party; // set owner, set active
  bullet[currentBullet].pos = pos, bullet[currentBullet].speed = speed;

  if (++currentBullet == MAXBULLET) currentBullet = 0;
}

// Tank::Tick - update single tank
void Tank::Tick()
{
  if (!(flags & ACTIVE)) // dead tank
  {
    smoke.xpos = (int)pos.x, smoke.ypos = (int)pos.y;
    return smoke.Tick();
  }

  // Complexity: O(n*k) (tanks*mountains)
  float2 force = Normalize(target - pos);
  // evade mountain peaks
  for (unsigned int i = 0; i < 16; i++)
  {
    float2 d(pos.x - peakx[i], pos.y - peaky[i]);
    float sd = (d.x * d.x + d.y * d.y) * 0.2f;
    if (sd < 1500) force += d * 0.03f * (peakh[i] / sd);
  }

  // evade other tanks
  quadTree->Collide(this, force);

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
  speed += force, speed = Normalize(speed), pos += speed * maxspeed * 0.5f;

  // shoot, if reloading completed
  if (--reloading >= 0) return;

  unsigned int start = 0, end = MAXP1;
  if (flags & P1) start = MAXP1, end = MAXP1 + MAXP2;

  for (unsigned int i = start; i < end; i++) if (game->m_Tank[i]->flags & ACTIVE)
  {
    float2 d = game->m_Tank[i]->pos - pos;
    if ((Length(d) < 100) && (Dot(Normalize(d), speed) > 0.99999f))
    {
      Fire(flags & (P1 | P2), pos, speed); // shoot
      reloading = 200; // and wait before next shot is ready
      break;
    }
  }
}

// Game::Init - Load data, setup playfield
void Game::Init()
{
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
  m_Tank = new Tank*[MAXP1 + MAXP2];
  m_P1Sprite = new Sprite(new Surface("testdata/p1tank.tga"), 1, Sprite::FLARE);
  m_P2Sprite = new Sprite(new Surface("testdata/p2tank.tga"), 1, Sprite::FLARE);
  m_PXSprite = new Sprite(new Surface("testdata/deadtank.tga"), 1, Sprite::BLACKFLARE);
  m_Smoke = new Sprite(new Surface("testdata/smoke.tga"), 10, Sprite::FLARE);
  // create blue tanks
  for (unsigned int i = 0; i < MAXP1; i++)
  {
    Tank* t = m_Tank[i] = new Tank();
    t->pos = float2((float)((i % 5) * 20), (float)((i / 5) * 20 + 50));
    t->target = float2(SCRWIDTH, SCRHEIGHT); // initially move to bottom right corner
    t->speed = float2(0, 0), t->flags = Tank::ACTIVE | Tank::P1, t->maxspeed = (i < (MAXP1 / 2)) ? 0.65f : 0.45f;
  }
  // create red tanks
  for (unsigned int i = 0; i < MAXP2; i++)
  {
    Tank* t = m_Tank[i + MAXP1] = new Tank();
    t->pos = float2((float)((i % 12) * 20 + 900), (float)((i / 12) * 20 + 600));
    t->target = float2(424, 336); // move to player base
    t->speed = float2(0, 0), t->flags = Tank::ACTIVE | Tank::P2, t->maxspeed = 0.3f;
  }
  game = this; // for global reference
  m_LButton = m_PrevButton = false;
}

// Game::DrawTanks - draw the tanks
void Game::DrawTanks()
{
  TimerRDTSC timer;

  timer.Start();
  for (unsigned int i = 0; i < MAXP1; i++)
  {
    if (!(m_Tank[i]->flags & Tank::ACTIVE)) continue;

    int xMin = m_Tank[i]->pos.x - 20;
    int yMin = m_Tank[i]->pos.y - 20;
    int xMax = xMin + 40;
    int yMax = yMin + 40;
    xMin = MAX(xMin, 0);
    yMin = MAX(yMin, 0);
    xMax = MIN(xMax, SCRWIDTH);
    yMax = MIN(yMax, SCRHEIGHT);

    for (int y = yMin; y < yMax; y++)
    {
      for (int x = xMin; x < xMax; x++)
      {
        float dx = m_Tank[i]->pos.x - x;
        float dy = m_Tank[i]->pos.y - y;
        float dist = dx * dx + dy * dy;

        if (dist < 400.f) // Sqdst of 20
          m_Surface->GetBuffer()[x + y * SCRWIDTH] |= 0x660000;
      }
    }
  }
  timer.Stop();
  auto timing = timer.Interval();

  char glowTimingsStr[128];
  sprintf(glowTimingsStr, "- Glow: %03i", timing);
  m_Surface->Print(glowTimingsStr, 20, 70, 0xffff00);

  timer.Start();
  for (unsigned int i = 0; i < (MAXP1 + MAXP2); i++)
  {
    Tank* t = m_Tank[i];
    float x = t->pos.x, y = t->pos.y;
    float2 p1(x + 70 * t->speed.x + 22 * t->speed.y, y + 70 * t->speed.y - 22 * t->speed.x);
    float2 p2(x + 70 * t->speed.x - 22 * t->speed.y, y + 70 * t->speed.y + 22 * t->speed.x);
    if (!(m_Tank[i]->flags & Tank::ACTIVE))
      m_PXSprite->Draw((int)x - 4, (int)y - 4, m_Surface); // draw dead tank
    else if (t->flags & Tank::P1) // draw blue tank
    {
      m_P1Sprite->Draw((int)x - 4, (int)y - 4, m_Surface);
      m_Surface->Line(x, y, x + 8 * t->speed.x, y + 8 * t->speed.y, 0x4444ff);
    }
    else // draw red tank
    {
      m_P2Sprite->Draw((int)x - 4, (int)y - 4, m_Surface);
      m_Surface->Line(x, y, x + 8 * t->speed.x, y + 8 * t->speed.y, 0xff4444);
    }

    // Drawing tracks if in screen. 
    if ((x >= 0) && (x < SCRWIDTH) && (y >= 0) && (y < SCRHEIGHT))
      m_Backdrop->GetBuffer()[(int)x + (int)y * SCRWIDTH] =
      SubBlend(m_Backdrop->GetBuffer()[(int)x + (int)y * SCRWIDTH], 0x030303); // tracks
  }
  timer.Stop();
  timing = timer.Interval();
  char drawTimingsStr[128];
  sprintf(drawTimingsStr, "- The rest: %03i", timing);
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
    for (unsigned int i = 0; i < MAXP1; i++) m_Tank[i]->target = float2((float)m_MouseX, (float)m_MouseY);
    m_Surface->Line(0, (float)m_MouseY, SCRWIDTH - 1, (float)m_MouseY, 0xffffff);
    m_Surface->Line((float)m_MouseX, 0, (float)m_MouseX, SCRHEIGHT - 1, 0xffffff);
  }
  m_PrevButton = m_LButton;
}

// Game::Tick - main game loop
void Game::Tick(float a_DT)
{
  POINT p;
  GetCursorPos(&p);
  ScreenToClient(FindWindow(NULL, "Template"), &p);
  m_LButton = (GetAsyncKeyState(VK_LBUTTON) != 0), m_MouseX = p.x, m_MouseY = p.y;

  TimerRDTSC timer;

  // Drawing backdrop
  m_Backdrop->CopyTo(m_Surface, 0, 0);

  // Update tanks
  timer.Start();
  quadTree = new QuadTree(0, SCRWIDTH, 0, SCRHEIGHT, 0);
  for (unsigned int i = 0; i < (MAXP1 + MAXP2); i++)
    quadTree->Add(m_Tank[i]);
  for (unsigned int i = 0; i < (MAXP1 + MAXP2); i++)
    m_Tank[i]->Tick();
  delete quadTree;
  timer.Stop();
  auto utTiming = timer.Interval();

  // Update bullets
  timer.Start();
  for (unsigned int i = 0; i < MAXBULLET; i++)
    bullet[i].Tick();
  timer.Stop();
  auto blTiming = timer.Interval();

  // Draw tanks
  timer.Start();
  DrawTanks();
  timer.Stop();
  auto dtTiming = timer.Interval();

  PlayerInput();

  //Timings:
  char utTimingsStr[128];
  sprintf(utTimingsStr, "Update Tanks: %03i", utTiming);
  m_Surface->Print(utTimingsStr, 10, 40, 0xffff00);

  char blTimingsStr[128];
  sprintf(blTimingsStr, "Update Bullets: %03i", blTiming);
  m_Surface->Print(blTimingsStr, 10, 50, 0xffff00);

  char dtTimingsStr[128];
  sprintf(dtTimingsStr, "Draw Tanks: %03i", dtTiming);
  m_Surface->Print(dtTimingsStr, 10, 60, 0xffff00);

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