// Template, major revision 5
// IGAD/NHTV - Jacco Bikker - 2006-2013

#pragma once

namespace Tmpl8 {

#include "emmintrin.h"
#include <smmintrin.h>

typedef unsigned long Pixel;
const unsigned int REDMASK = (0xff0000);
const unsigned int GREENMASK = (0x00ff00);
const unsigned int BLUEMASK = (0x0000ff);

static const __m128i REDMASK4 = _mm_set1_epi32(REDMASK);
static const __m128i GREENMASK4 = _mm_set1_epi32(GREENMASK);
static const __m128i BLUEMASK4 = _mm_set1_epi32(BLUEMASK);

static __m128i AddBlend4(const __m128i &a_Color1, const __m128i &a_Color2)
{

  //const unsigned int r = (a_Color1 & REDMASK) + (a_Color2 & REDMASK);
  //const unsigned int g = (a_Color1 & GREENMASK) + (a_Color2 & GREENMASK);
  //const unsigned int b = (a_Color1 & BLUEMASK) + (a_Color2 & BLUEMASK);
  const __m128i r = _mm_add_epi32(_mm_and_si128(a_Color1, REDMASK4), _mm_and_si128(a_Color2, REDMASK4));
  const __m128i g = _mm_add_epi32(_mm_and_si128(a_Color1, GREENMASK4), _mm_and_si128(a_Color2, GREENMASK4));
  const __m128i b = _mm_add_epi32(_mm_and_si128(a_Color1, BLUEMASK4), _mm_and_si128(a_Color2, BLUEMASK4));

  //const unsigned r1 = (r & REDMASK) | (REDMASK * (r >> 24));
  //const unsigned g1 = (g & GREENMASK) | (GREENMASK * (g >> 16));
  //const unsigned b1 = (b & BLUEMASK) | (BLUEMASK * (b >> 8));

  const __m128i r1 = _mm_or_si128(_mm_and_si128(r, REDMASK4), _mm_mullo_epi32(REDMASK4, _mm_srai_epi32(r, 24)));
  const __m128i g1 = _mm_or_si128(_mm_and_si128(g, GREENMASK4), _mm_mullo_epi32(GREENMASK4, _mm_srai_epi32(g, 16)));
  const __m128i b1 = _mm_or_si128(_mm_and_si128(b, BLUEMASK4), _mm_mullo_epi32(BLUEMASK4, _mm_srai_epi32(b, 8)));

  //return (r1 + g1 + b1);
  return _mm_add_epi32(_mm_add_epi32(r1, g1), b1);
}
/*
static __m128 AddBlend4(__m128 a_Color1, __m128 a_Color2)
{
  __m128 col1Redmasked = _mm_and_ps(a_Color1, REDMASK4), col2Redmasked = _mm_and_ps(a_Color2, REDMASK4);
  const __m128 r = *(__m128*)&_mm_add_epi32(*(__m128i*)&col1Redmasked, *(__m128i*)&col2Redmasked);

  __m128 col1Greenmasked = _mm_and_ps(a_Color1, GREENMASK4), col2Greenmasked = _mm_and_ps(a_Color2, GREENMASK4);
  const __m128 g = *(__m128*)&_mm_add_epi32(*(__m128i*)&col1Greenmasked, *(__m128i*)&col2Greenmasked);

  __m128 col1Bluemasked = _mm_and_ps(a_Color1, BLUEMASK4), col2Bluemasked = _mm_and_ps(a_Color2, BLUEMASK4);
  const __m128 b = *(__m128*)&_mm_add_epi32(*(__m128i*)&col1Bluemasked, *(__m128i*)&col2Bluemasked);

  const __m128 r1 = _mm_or_ps(_mm_and_ps(r, REDMASK4), *(__m128*)&_mm_mul_epu32(*(__m128i*)&REDMASK4, _mm_srli_epi32(*(__m128i*)&REDMASK4, 24))); 
  const __m128 g1 = _mm_or_ps(_mm_and_ps(g, GREENMASK4), *(__m128*)&_mm_mul_epu32(*(__m128i*)&GREENMASK4, _mm_srli_epi32(*(__m128i*)&GREENMASK4, 16)));
  const __m128 b1 = _mm_or_ps(_mm_and_ps(b, BLUEMASK4), *(__m128*)&_mm_mul_epu32(*(__m128i*)&BLUEMASK4, _mm_srli_epi32(*(__m128i*)&BLUEMASK4, 8)));
  
  return _mm_add_ps(_mm_add_ps(r1, g1), b1);
}*/

inline Pixel AddBlend(Pixel a_Color1, Pixel a_Color2)
{
  const unsigned int r = (a_Color1 & REDMASK) + (a_Color2 & REDMASK);
  const unsigned int g = (a_Color1 & GREENMASK) + (a_Color2 & GREENMASK);
  const unsigned int b = (a_Color1 & BLUEMASK) + (a_Color2 & BLUEMASK);
  const unsigned r1 = (r & REDMASK) | (REDMASK * (r >> 24));
  const unsigned g1 = (g & GREENMASK) | (GREENMASK * (g >> 16));
  const unsigned b1 = (b & BLUEMASK) | (BLUEMASK * (b >> 8));
  return (r1 + g1 + b1);
}

// subtractive blending
inline Pixel SubBlend( Pixel a_Color1, Pixel a_Color2 )
{
	int red = (a_Color1 & REDMASK) - (a_Color2 & REDMASK);
	int green = (a_Color1 & GREENMASK) - (a_Color2 & GREENMASK);
	int blue = (a_Color1 & BLUEMASK) - (a_Color2 & BLUEMASK);
	if (red < 0) red = 0;
	if (green < 0) green = 0;
	if (blue < 0) blue = 0;
	return (Pixel)(red + green + blue);
}

inline Pixel ScaleColor( Pixel a_Color, unsigned int a_Scale )
{
	unsigned int rb = (((a_Color & (REDMASK|BLUEMASK)) * a_Scale) >> 8) & (REDMASK|BLUEMASK);
	unsigned int g = (((a_Color & GREENMASK) * a_Scale) >> 8) & GREENMASK;
	return rb + g;
}

class Color
{
public:
	Color() : r( 0.0f ), g( 0.0f ), b( 0.0f ), a( 0.0f ) {};
	Color( float a_R, float a_G, float a_B ) : r( a_R ), g( a_G ), b( a_B ), a( 0.0f ) {};
	void Set( float a_R, float a_G, float a_B ) { r = a_R; g = a_G; b = a_B; a = 0; }
	void operator += ( const Color& a_V ) { r += a_V.r; g += a_V.g; b += a_V.b; }
	void operator += ( Color* a_V ) { r += a_V->r; g += a_V->g; b += a_V->b; }
	void operator -= ( const Color& a_V ) { r -= a_V.r; g -= a_V.g; b -= a_V.b; }
	void operator -= ( Color* a_V ) { r -= a_V->r; g -= a_V->g; b -= a_V->b; }
	void operator *= ( const float f ) { r *= f; g *= f; b *= f; }
	void operator *= ( const Color& a_V ) { r *= a_V.r; g *= a_V.g; b *= a_V.b; }
	void operator *= ( Color* a_V ) { r *= a_V->r; g *= a_V->g; b *= a_V->b; }
	Color operator- () const { return Color( -r, -g, -b ); }
	friend Color operator + ( const Color& v1, const Color& v2 ) { return Color( v1.r + v2.r, v1.g + v2.g, v1.b + v2.b ); }
	friend Color operator - ( const Color& v1, const Color& v2 ) { return Color( v1.r - v2.r, v1.g - v2.g, v1.b - v2.b ); }
	friend Color operator + ( const Color& v1, Color* v2 ) { return Color( v1.r + v2->r, v1.g + v2->g, v1.b + v2->b ); }
	friend Color operator - ( const Color& v1, Color* v2 ) { return Color( v1.r - v2->r, v1.g - v2->g, v1.b - v2->b ); }
	friend Color operator * ( const Color& v, const float f ) { return Color( v.r * f, v.g * f, v.b * f ); }
	friend Color operator * ( const Color& v1, const Color& v2 ) { return Color( v1.r * v2.r, v1.g * v2.g, v1.b * v2.b ); }
	friend Color operator * ( const float f, const Color& v ) { return Color( v.r * f, v.g * f, v.b * f ); }
	friend Color operator / ( const Color& v, const float f ) { float r = 1.0f / f; return Color( v.r * r, v.g * r, v.b * r ); }
	union
	{
		struct { float b, g, r, a; };
		float cell[4];
		__m128 rgba;
		__m128i irgba;
	};
};

class Surface
{
	enum
	{
		OWNER = 1
	};

public:
	// constructor / destructor
	Surface( int a_Width, int a_Height, Pixel* a_Buffer, int a_Pitch );
	Surface( int a_Width, int a_Height );
	Surface( char* a_File );
	~Surface();
	// member data access
	Pixel* GetBuffer() { return m_Buffer; }
	void SetBuffer( Pixel* a_Buffer ) { m_Buffer = a_Buffer; }
	int GetWidth() { return m_Width; }
	int GetHeight() { return m_Height; }
	int GetPitch() { return m_Pitch; }
	void SetPitch( int a_Pitch ) { m_Pitch = a_Pitch; }
	void ParseHeader( unsigned char* a_Header );
	// Special operations
	void InitCharset();
	void SetChar( int c, char* c1, char* c2, char* c3, char* c4, char* c5 );
	void Centre( char* a_String, int y1, Pixel color );
	void Print( char* a_String, int x1, int y1, Pixel color );
	void Clear( Pixel a_Color );
	void Line( float x1, float y1, float x2, float y2, Pixel color );
	void AddLine( float x1, float y1, float x2, float y2, Pixel color );
	void ThickLine( int ax1, int ay1, int ax2, int ay2, Pixel c );
	void Plot( int x, int y, Pixel c );
	void LoadImage( char* a_File );
	void CopyTo( Surface* a_Dst, int a_X, int a_Y );
	void BlendCopyTo( Surface* a_Dst, int a_X, int a_Y );
	void ScaleColor( unsigned int a_Scale );
	void Box( int x1, int y1, int x2, int y2, Pixel color );
	void Bar( int x1, int y1, int x2, int y2, Pixel color );
  void Resize(Surface* a_Orig);
  int m_Width, m_Height, m_Pitch;
private:
	// Attributes
	Pixel* m_Buffer;	
	// Static attributes for the buildin font
	static char s_Font[51][5][5];	
	static int s_Transl[256];		
};

class Sprite
{
public:
	// Sprite flags
	enum
	{
		FLARE		= (1<< 0),
		OPFLARE		= (1<< 1),	
		FLASH		= (1<< 4),	
		DISABLED	= (1<< 6),	
		GMUL		= (1<< 7),
		BLACKFLARE	= (1<< 8),	
		BRIGHTEST   = (1<< 9),
		RFLARE		= (1<<12),
		GFLARE		= (1<<13),
		NOCLIP		= (1<<14)
	};
	
	// Structors
	Sprite( Surface* a_Surface, unsigned int a_NumFrames );
	Sprite( Surface* a_Surface, unsigned int a_NumFrames, unsigned int a_Flags );
	~Sprite();
	// Methods
	void Draw( int a_X, int a_Y, Surface* a_Target = 0 );
	void DrawScaled( int a_X, int a_Y, int a_Width, int a_Height, Surface* a_Target );
	void SetFlags( unsigned int a_Flags ) { m_Flags = a_Flags; }
	void SetFrame( unsigned int a_Index ) { m_CurrentFrame = a_Index; }
	unsigned int GetFlags() const { return m_Flags; }
	int GetWidth() { return m_Width; }
	int GetHeight() { return m_Height; }
	Pixel* GetBuffer() { return m_Surface->GetBuffer(); }	
	unsigned int Frames() { return m_NumFrames; }
  Surface* GetSurface() { return m_Surface; }
private:
	// Methods
	void InitializeStartData();
	// Attributes
	int m_Width, m_Height, m_Pitch;
	unsigned int m_NumFrames;          
	unsigned int m_CurrentFrame;       
	unsigned int m_Flags;
	unsigned int** m_Start;
	Surface* m_Surface;
  Pixel** m_FramePtrs;
};

class Font
{
public:
	Font();
	Font( char* a_File, char* a_Chars );
	~Font();
	void Print( Surface* a_Target, char* a_Text, int a_X, int a_Y, bool clip = false );
	void Centre( Surface* a_Target, char* a_Text, int a_Y );
	int Width( char* a_Text );
	int Height() { return m_Surface->GetHeight(); }
	void YClip( int y1, int y2 ) { m_CY1 = y1; m_CY2 = y2; }
private:
	Surface* m_Surface;
	int* m_Offset, *m_Width, *m_Trans, m_Height, m_CY1, m_CY2;
};

class AlignedSprite
{
public:
  AlignedSprite(Surface* a_Surface);

  __declspec(align(16)) union{ __m128i alignmentBuffers4[16][12][3]; Pixel alignmentBuffers[16][12][12]; };
  void FillBuffer(Pixel* a_buffer, int a_x, int a_y);
  
  void Draw(int a_x, int a_y, Surface* a_target);
};

}; // namespace Tmpl8