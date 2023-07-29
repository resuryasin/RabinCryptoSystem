#include <iostream>
#include <algorithm>
#include <iterator>
#include <string>
#include <cstring>
#include <cmath>
#include "mini-gmp.h"
#include "mini-gmp.c"
#include "variables.cpp"

using namespace std;

void print_vec(vector < int > v) {
  copy(v.begin(), v.end(), ostream_iterator < int > (cout, " "));
}

vector < int > dec2bin(long m) {
  vector < int > l;
  long temp;
  while (m > 0) {
    temp = m % 2;
    l.insert(l.begin(), temp);
    m = m / 2;
  }
  return l;
}

long bin2dec(vector < int > m) {
  long power, decimal = 0;
  for (int i = m.size() - 1; i >= 0; i--) {
    power = m.size() - i - 1;
    decimal += m[i] * pow(2, power);
  }
  return decimal;
}

vector < int > doublebin(vector < int > m) {
  vector < int > l;
  for (int i = 0; i < m.size(); i++) {
    l.push_back(m[i]);
  }
  for (int i = 0; i < m.size(); i++) {
    l.push_back(m[i]);
  }
  return l;
}

long kongruen(long m, long n) {
  long k;
  k = (m - (m % n)) / n;
  return k;
}

long encrypt(long m, long n) {
  long c;
  mpz_t tmp;
  mpz_init(tmp);
  mpz_set_ui(tmp, m);
  mpz_pow_ui(tmp, tmp, 2);
  mpz_mod_ui(tmp, tmp, n);
  c = mpz_get_ui(tmp);
  mpz_clear(tmp);
  return c;
}

long2 euclidex(long p, long q) {
  long fx, fy, tempp, tempq;
  tempp = p;
  tempq = q;
  long x = 0, y = 1;
  long lastx = 1, lasty = 0, temp;
  while (tempq != 0) {
    long z1 = tempp / tempq;
    long z2 = tempp % tempq;
    tempp = tempq;
    tempq = z2;

    temp = x;
    x = lastx - z1 * x;
    lastx = temp;

    temp = y;
    y = lasty - z1 * y;
    lasty = temp;
  }
  fx = lastx;
  fy = lasty;
  long2 euxy = {
    fx,
    fy
  };
  return euxy;
}

long2 quampq(long c, long p, long q) {
  long2 mpq, pq;
  pq.x = (p + 1) / 4;
  pq.y = (q + 1) / 4;
  mpz_t mp, mp1, mp2;
  mpz_init(mp);
  mpz_init(mp1);
  mpz_init(mp2);
  mpz_set_ui(mp, c);
  mpz_pow_ui(mp1, mp, pq.x);
  mpz_pow_ui(mp2, mp, pq.y);
  mpz_mod_ui(mp1, mp1, p);
  mpz_mod_ui(mp2, mp2, q);
  mpq.x = mpz_get_ui(mp1);
  mpq.y = mpz_get_ui(mp2);
  mpz_clear(mp);
  mpz_clear(mp1);
  mpz_clear(mp2);
  return mpq;
}

long4 crt(long2 eu, long2 mpq, long p, long q, long n) {
  long4 cr;
  cr.w = (eu.x * p * mpq.y + eu.y * q * mpq.x) % n;
  cr.x = (-cr.w) % n;
  cr.y = (eu.y * q * mpq.x - eu.x * p * mpq.y) % n;
  cr.z = (-cr.y) % n;
  if (cr.w < 0) cr.w = n + cr.w;
  if (cr.x < 0) cr.x = n + cr.x;
  if (cr.y < 0) cr.y = n + cr.y;
  if (cr.z < 0) cr.z = n + cr.z;
  return cr;
}

long4 krstu(long k, long4 rstu, long n) {
  long4 l;
  l.w = (k * n) + rstu.w;
  l.x = (k * n) + rstu.x;
  l.y = (k * n) + rstu.y;
  l.z = (k * n) + rstu.z;
  return l;
}

long decrypt(vector < int > c) {
  vector < int > l1, l2;
  long result = 0, temp;
  if (c.size() % 2 == 0) {
    temp = c.size() / 2;
  } else {
    temp = (c.size() + 1) / 2;
  }
  for (int i = 0; i < temp; i++) {
    l1.push_back(c[i]);
  }
  for (int j = c.size() - 1; j >= temp; j--) {
    l2.insert(l2.begin(), c[j]);
  }
  int k;
  for (k = 0; k < temp; k++) {
    if (l1[k] != l2[k]) {
      break;
    }
  }
  if (k == temp) {
    result = bin2dec(l1);
  }
  return result;
}

long new_message(long m) {
  long new_m;
  vector < int > m_bin, m_double_bin;
  m_bin = dec2bin(m);
  m_double_bin = doublebin(m_bin);
  new_m = bin2dec(m_double_bin);
  return new_m;
}

int pdekripsi(long c, long p, long q, long k, long n) {
  long cma;
  long2 eu, mpq;
  long4 rstu, nrstu, cm;
  vector4 cbin;

  eu = euclidex(p, q);
  mpq = quampq(c, p, q);
  rstu = crt(eu, mpq, p, q, n);
  nrstu = krstu(k, rstu, n);
  cbin.a = dec2bin(nrstu.w);
  cbin.b = dec2bin(nrstu.x);
  cbin.c = dec2bin(nrstu.y);
  cbin.d = dec2bin(nrstu.z);

  cm.w = decrypt(cbin.a);
  cm.x = decrypt(cbin.b);
  cm.y = decrypt(cbin.c);
  cm.z = decrypt(cbin.d);
  if (cm.w != 0) {
    cma = cm.w;
  } else if (cm.x != 0) {
    cma = cm.x;
  } else if (cm.y != 0) {
    cma = cm.y;
  } else if (cm.z != 0) {
    cma = cm.z;
  } else {
    cout << "Error" << endl;
  }
  return cma;
}
int main() {
  long p, q, n, k, c, cma;
  vector < int > mbin, mbin2;
  long2 eu, mpq, test;
  long4 rstu, nrstu, cm;
  vector4 cbin;

  p = 23;
  q = 11;
  n = p * q;

  string msg = "Keamanan Data";
  int len = strlen(msg.c_str());
  vector < int > vmsg, vmsgbaru, vk, vc, vpt;
  for (int i = 0; i < len; i++) {
    vmsg.push_back((msg[i]));
  }
  for (int i = 0; i < len; i++) {
    vmsgbaru.push_back(new_message(msg[i]));
  }
  for (int i = 0; i < len; i++) {
    vk.push_back(kongruen(vmsgbaru[i], n));
  }
  for (int i = 0; i < len; i++) {
    vc.push_back(encrypt(vmsgbaru[i], n));
  }
  for (int i = 0; i < len; i++) {
    vpt.push_back(pdekripsi(vc[i], p, q, vk[i], n));
  }

  cout << "Pesan Asli : " << msg << endl << "Kunci Publik: " << n << endl << endl << "Pesan Asli (ASCII): ";
  print_vec(vmsg);
  cout << endl << "Pesan Diduplikasi: ";
  print_vec(vmsgbaru);
  cout << endl << "Pesan Dienkripsi: ";
  print_vec(vc);
  cout << endl << "Chiper Text: ";
  copy(vc.begin(), vc.end(), ostream_iterator < char > (cout, ""));
  cout << endl << endl << "Kongruen: ";
  print_vec(vk);
  cout << endl << "Pesan Didekripsi: ";
  print_vec(vpt);
  cout << endl << "Plain Text: ";
  copy(vpt.begin(), vpt.end(), ostream_iterator < char > (cout, ""));
}