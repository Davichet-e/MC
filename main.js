// Bilbiography:
// https://codereview.stackexchange.com/questions/124225/2d-colliding-disks-in-javascript-ideal-billiard-balls

const eps = 1e-6;

/**
 * Vec2
 *
 * Clase de ayuda para operaciones con vectores
 */
class Vec2 {
  constructor(x, y) {
    this.x = x;
    this.y = y;
    this.c = 2;
  }
  plus(p) {
    return new Vec2(this.x + p.x, this.y + p.y);
  }
  minus(p) {
    return new Vec2(this.x - p.x, this.y - p.y);
  }
  scale(v) {
    return new Vec2(this.x * v, this.y * v);
  }
  times(v) {
    return new Vec2(this.x * v, this.y * v);
  }
  length() {
    return Math.sqrt(this.x * this.x + this.y * this.y);
  }
  dot(other) {
    return this.x * other.x + this.y * other.y;
  }
  transform(t) {
    return new Vec2(
      t.m[0] * this.x + t.m[1] * this.y + t.m[2],
      t.m[3] * this.x + t.m[4] * this.y + t.m[5]
    );
  }
  distance(other) {
    return this.minus(other).length();
  }
  unit() {
    return this.scale(1 / (this.length() || 1));
  }
  neg() {
    return this.scale(-1);
  }
  project(other) {
    other = other.unit();
    return other.scale(this.dot(other));
  }
  reflect(other) {
    const projection = this.project(other);
    const residual = this.minus(projection);
    return this.minus(residual.times(2));
  }
}

/**
 * Vec2
 *
 * Clase de ayuda para operaciones con elipses
 */
class Ellipse {
  constructor(a, b, c, d, e, f) {
    if (a === undefined) {
      this.a = 1;
      this.b = 1;
      this.c = 0;
      this.d = 0;
      this.e = 0;
      this.f = -1;
    } else {
      this.a = a;
      this.b = b;
      this.c = c;
      this.d = d;
      this.e = e;
      this.f = f;
    }
    this.gradient = new Transform(2 * a, c, d, c, 2 * b, e);
    this.center = this.gradient.inverse().apply(new Vec2(0, 0));
  }
  transform(t) {
    const i = t.inverse();
    const [m00, m01, m02, m10, m11, m12] = i.m;
    const aa = this.a * m00 * m00 + this.b * m10 * m10 + this.c * m00 * m10;
    const bb = this.a * m01 * m01 + this.b * m11 * m11 + this.c * m01 * m11;
    const cc =
      2 * this.a * m00 * m01 +
      2 * this.b * m10 * m11 +
      this.c * (m00 * m11 + m01 * m10);
    const dd =
      2 * this.a * m00 * m02 +
      2 * this.b * m10 * m12 +
      this.c * (m00 * m12 + m02 * m10) +
      this.d * m00 +
      this.e * m10;
    const ee =
      2 * this.a * m01 * m02 +
      2 * this.b * m11 * m12 +
      this.c * (m01 * m12 + m02 * m11) +
      this.d * m01 +
      this.e * m11;
    const ff =
      this.a * m02 * m02 +
      this.b * m12 * m12 +
      this.c * m02 * m12 +
      this.d * m02 +
      this.e * m12 +
      this.f;
    return new Ellipse(aa, bb, cc, dd, ee, ff);
  }
  lineIntersection(c, p) {
    const pc = p.minus(c);
    const u2 =
      this.a * Math.pow(pc.x, 2) +
      this.b * Math.pow(pc.y, 2) +
      this.c * pc.x * pc.y;
    const u1 =
      2 * this.a * c.x * pc.x +
      2 * this.b * c.y * pc.y +
      this.c * c.y * pc.x +
      this.c * c.x * pc.y +
      this.d * pc.x +
      this.e * pc.y;
    const u0 =
      this.a * Math.pow(c.x, 2) +
      this.b * Math.pow(c.y, 2) +
      this.c * c.x * c.y +
      this.d * c.x +
      this.e * c.y +
      this.f;
    const result = [];
    const sols = quadratic(u2, u1, u0);
    if (isNaN(sols.root1) || isNaN(sols.root2)) {
      return result;
    }
    result.push(c.plus(pc.scale(sols.root1)));
    result.push(c.plus(pc.scale(sols.root2)));
    return result;
  }
}

function withinEps(v) {
  return Math.abs(v) < eps;
}

function quadratic(a, b, c) {
  if (a === 0) {
    return { root1: -c / b, root2: -c / b };
  }
  if (b === 0 && c === 0) {
    return { root1: 0, root2: 0 };
  }
  let discriminant = b * b - 4 * a * c;
  if (discriminant < 0 && withinEps(discriminant)) {
    discriminant = 0;
  }
  if (discriminant === 0) {
    return { root1: -b / (2 * a), root2: -b / (2 * a) };
  }
  const d = Math.sqrt(discriminant);

  if (b >= 0) {
    return { root1: (-b - d) / (2 * a), root2: (2 * c) / (-b - d) };
  } else {
    return { root1: (2 * c) / (-b + d), root2: (-b + d) / (2 * a) };
  }
}

class Transform {
  constructor(m11, m12, tx, m21, m22, ty) {
    if (m12 === undefined) {
      this.m = m11;
    } else {
      this.m = new Float32Array([m11, m12, tx, m21, m22, ty, 0, 0, 1]);
    }
  }
  static rotate(theta) {
    const s = Math.sin(theta);
    const c = Math.cos(theta);
    return new Transform(c, -s, 0, s, c, 0);
  }
  static translate(tx, ty) {
    if (ty === undefined) {
      return new Transform(1, 0, tx.x, 0, 1, tx.y);
    }
    return new Transform(1, 0, tx, 0, 1, ty);
  }
  static scale(sx, sy) {
    return new Transform(sx, 0, 0, 0, sy, 0);
  }
  det() {
    return this.m[0] * this.m[4] - this.m[1] * this.m[3];
  }
  apply(other) {
    return other.transform(this);
  }
  get(i, j) {
    return this.m[i * 3 + j];
  }
  set(i, j, v) {
    this.m[i * 3 + j] = v;
  }
  inverse() {
    const d = 1.0 / this.det();
    const t = new Transform(
      d * this.get(1, 1),
      -d * this.get(0, 1),
      0,
      -d * this.get(1, 0),
      d * this.get(0, 0),
      0
    );
    const v = t.apply(new Vec2(this.get(0, 2), this.get(1, 2)));
    t.set(0, 2, -v.x);
    t.set(1, 2, -v.y);
    return t;
  }
}

function interseccionLineas(p1, p2, p3, p4) {
  function det2(a, b, c, d) {
    return a * d - b * c;
  }
  const x12 = p1.x - p2.x;
  const x34 = p3.x - p4.x;
  const y12 = p1.y - p2.y;
  const y34 = p3.y - p4.y;
  const t1 = det2(p1.x, p1.y, p2.x, p2.y);
  const t2 = det2(p3.x, p3.y, p4.x, p4.y);
  const t3 = det2(x12, y12, x34, y34);
  return new Vec2(det2(t1, x12, t2, x34) / t3, det2(t1, y12, t2, y34) / t3);
}

function hazBola(origin, velocity, originT) {
  return {
    originT: originT ?? 0,
    origin,
    velocity,
  };
}

function botaBola(shape, ball) {
  const intersection = shape.intersect(ball);
  if (
    intersection.t === Infinity ||
    isNaN(intersection.t) ||
    intersection.t <= ball.originT
  )
    return undefined;

  const normal = shape.normalAt(intersection.point);
  if (normal.dot(ball.velocity.unit()) >= 0) return undefined;

  return hazBola(
    intersection.point,
    ball.velocity.reflect(normal).neg(),
    intersection.t
  );
}

function medioPlano(point, normal) {
  normal = normal.unit();
  const tangent = new Vec2(normal.y, -normal.x);
  const p1 = point.plus(tangent),
    p2 = point.minus(tangent);

  return {
    intersect(ball) {
      const result = interseccionLineas(
        p1,
        p2,
        ball.origin,
        ball.origin.plus(ball.velocity)
      );
      let distT = result.distance(ball.origin) / ball.velocity.length();
      if (normal.dot(ball.velocity) > 0) distT = -distT;
      let intersectT = ball.originT + distT;
      if (intersectT < ball.originT) {
        intersectT = Infinity;
      }
      return {
        t: intersectT,
        point: result,
      };
    },
    normalAt() {
      return normal;
    },
  };
}

function medioCirculo(center, radius, normal) {
  const ellipse = new Ellipse()
    .transform(Transform.scale(radius, radius))
    .transform(Transform.translate(center.x, center.y));

  normal = normal.unit();
  return {
    intersect(ball) {
      const is = ellipse.lineIntersection(
        ball.origin,
        ball.origin.plus(ball.velocity)
      );
      if (is.length !== 2) return undefined;

      const d1 = is[0].minus(ball.origin).unit();
      const d2 = is[1].minus(ball.origin).unit();
      const ballDirection = ball.velocity.unit();
      let result;
      // the intersections are only valid if they belong to the half-circle
      if (
        is[0].minus(center).unit().dot(normal) < 0 &&
        d1.dot(this.normalAt(is[0])) < 0 &&
        d1.dot(ballDirection) > eps
      ) {
        result = is[0];
      } else if (
        is[1].minus(center).unit().dot(normal) < 0 &&
        d2.dot(this.normalAt(is[1])) < 0 &&
        d2.dot(ballDirection) > eps
      ) {
        result = is[1];
      } else return undefined;
      const distT = result.distance(ball.origin) / ball.velocity.length();
      const intersectT = ball.originT + distT;

      return {
        t: intersectT,
        point: result,
      };
    },
    normalAt(point) {
      return ellipse.gradient.apply(point).unit().neg();
    },
  };
}

function decideSiguienteBote(ball, escena) {
  const botes = [];
  escena.forEach((object) => {
    const interseccion = object.intersect(ball);
    if (interseccion === undefined) return;
    if (interseccion.t === Infinity) return;
    const bote = botaBola(object, ball);
    if (bote !== undefined) {
      botes.push(bote);
    }
  });
  botes.sort((bote1, bote2) => {
    if (bote1.originT < bote2.originT) return -1;
    else if (bote1.originT > bote2.originT) return 1;
    else return 0;
  });
  if (botes.length === 0) {
    console.error("Bug");
  }
  return botes[0];
}

const svg = d3
  .select("#contenedor-sensibilidad")
  .append("svg")
  .style("position", "absolute")
  .attr("width", 800)
  .attr("height", 400);

const escena = [
  medioPlano(new Vec2(0, 200), new Vec2(0, -1)),
  medioPlano(new Vec2(0, -200), new Vec2(0, 1)),
  medioCirculo(new Vec2(-200, 0), 200, new Vec2(1, 0)),
  medioCirculo(new Vec2(200, 0), 200, new Vec2(-1, 0)),
];

function añadeFondoBunimovich(sel) {
  sel
    .attr("transform", "translate(400, 200)")
    .append("path")
    .attr(
      "d",
      "M -200 -200 l 400 0 a 200 200 0 0 1 0 400 l -400 0 a 200 200 0 0 1 0 -400 "
    )
    .attr("fill", d3.lab(100, 400, 0));
}

const stadiumBg = svg.append("g");

d3.select("#contenedor-sensibilidad").style("position", "relative");

const canvas = d3
  .select("#contenedor-sensibilidad")
  .append("canvas")
  .style("position", "absolute")
  .attr("width", 800)
  .attr("height", 400)
  .style("opacity", 0.7)
  .node();

d3.select("#contenedor-sensibilidad")
  .append("div")
  .style("position", "relative")
  .style("height", "420px")
  .style("width", "800px");

const ctx = canvas.getContext("2d");

// Esto indica de donde nacen las bolas
const svgScene = svg.append("g").attr("transform", "translate(400, 200)");

let bolas = [];

function añadeBolas(n, centro, direccion, diferencia, velocidad) {
  const angleScale = d3
    .scaleLinear()
    .domain([0, n])
    .range([-diferencia, diferencia]);
  for (let i = 0; i < n; ++i) {
    const rot = Transform.rotate(angleScale(i + (0.5 * Math.random() - 0.25)));
    bolas.push({
      current: hazBola(centro, rot.apply(direccion).scale(velocidad)),
    });
  }
}

const canvasBola = svgScene.append("g");

document.getElementById("reinicio").addEventListener("click", () => {
  demuestraSensibilidad(escena, añadeFondoBunimovich);
});

let iteracion = 0;
const velocidad = 5;
function demuestraSensibilidad(escena, backgroundFunction) {
  ++iteracion;
  bolas = [];
  ctx.clearRect(0, 0, 800, 400);
  console.log(backgroundFunction);

  backgroundFunction(stadiumBg);

  añadeBolas(
    300,
    new Vec2(Math.random() * 100 - 50, Math.random() * 100 - 50),
    new Vec2(80, 80),
    0.00001,
    velocidad
  );

  const todasLasBolas = canvasBola.selectAll("*");
  todasLasBolas.transition();
  todasLasBolas.remove();

  const bolasVerdes = canvasBola
    .selectAll("*")
    .data(bolas)
    .enter()
    .append("circle")
    .each((ball) => {
      ball.next = decideSiguienteBote(ball.current, escena);
    })
    .attr("cx", function (ball) {
      return ball.current.origin.x;
    })
    .attr("cy", function (ball) {
      return ball.current.origin.y;
    })
    .attr("r", 5)
    .attr("fill", "green");
  bolasVerdes.call(siguienteTransicion(iteracion, escena, false));

  const bolasAzules = d3.select(bolasVerdes.nodes()[1]).attr("fill", "blue");
  bolasAzules.call(siguienteTransicion(iteracion, escena, true));
}

function siguienteTransicion(nDeIteracion, escena, pintaBolaBlanca) {
  return function (sel) {
    sel
      .transition()
      .ease(d3.easeLinear)
      .duration((bola) => {
        return (bola.next.originT - bola.current.originT) * 1000;
      })
      .attr("cx", (bola) => {
        return bola.next.origin.x;
      })
      .attr("cy", (bola) => {
        return bola.next.origin.y;
      })
      .on("end", function (bola) {
        try {
          ctx.beginPath();
          ctx.moveTo(bola.current.origin.x + 400, bola.current.origin.y + 200);
          ctx.lineTo(bola.next.origin.x + 400, bola.next.origin.y + 200);
          ctx.lineWidth = 3;
          ctx.strokeStyle = pintaBolaBlanca ? "white" : "black";
          ctx.stroke();
          bola.current = bola.next;
          bola.next = decideSiguienteBote(bola.current, escena);
          if (nDeIteracion == iteracion)
            d3.select(this).call(
              siguienteTransicion(nDeIteracion, escena, pintaBolaBlanca)
            );
        } catch (e) {
          console.error("");
        }
      });
  };
}

demuestraSensibilidad(escena, añadeFondoBunimovich);

///
const svg2 = d3
  .select("#contenedor-periodicidad")
  .append("svg")
  .style("position", "absolute")
  .attr("width", 800)
  .attr("height", 400);

const stadiumBg2 = svg2.append("g");

d3.select("#contenedor-periodicidad").style("position", "relative");

const canvas2 = d3
  .select("#contenedor-periodicidad")
  .append("canvas")
  .style("position", "absolute")
  .attr("width", 800)
  .attr("height", 400)
  .style("opacity", 0.7)
  .node();

d3.select("#contenedor-periodicidad")
  .append("div")
  .style("position", "relative")
  .style("height", "420px")
  .style("width", "800px");

const ctx2 = canvas2.getContext("2d");

// Esto indica de donde nacen las bolas
const svgScene2 = svg2.append("g").attr("transform", "translate(400, 200)");

let bolas2 = [];

function añadeBolas2(n, centro, direccion, diferencia, velocidad) {
  const angleScale = d3
    .scaleLinear()
    .domain([0, n])
    .range([-diferencia, diferencia]);
  for (let i = 0; i < n; ++i) {
    const rot = Transform.rotate(angleScale(i + (0.5 * Math.random() - 0.25)));
    bolas2.push({
      current: hazBola(centro, rot.apply(direccion).scale(velocidad)),
    });
  }
}

const canvasBola2 = svgScene2.append("g");

let iteracion2 = 0;
const velocidad2 = 5;
function demuestra2(escena, backgroundFunction) {
  ++iteracion2;
  bolas2 = [];
  ctx2.clearRect(0, 0, 800, 400);
  console.log(backgroundFunction);

  backgroundFunction(stadiumBg2);

  añadeBolas2(
    1,
    new Vec2(0, 0), // center
    new Vec2(0, -90), //direction
    0.00001,
    velocidad2
  );

  const todasLasBolas = canvasBola2.selectAll("*");
  todasLasBolas.transition();
  todasLasBolas.remove();

  const bolasVerdes = canvasBola2
    .selectAll("*")
    .data(bolas2)
    .enter()
    .append("circle")
    .each((ball) => {
      ball.next = decideSiguienteBote(ball.current, escena);
    })
    .attr("cx", function (ball) {
      return ball.current.origin.x;
    })
    .attr("cy", function (ball) {
      return ball.current.origin.y;
    })
    .attr("r", 15)
    .attr("fill", "green");
  bolasVerdes.call(siguienteTransicion2(iteracion2, escena, false));

  const bolasAzules = d3.select(bolasVerdes.nodes()[1]).attr("fill", "blue");
  bolasAzules.call(siguienteTransicion2(iteracion2, escena, true));
}

function siguienteTransicion2(nDeIteracion, escena, pintaBolaBlanca) {
  return function (sel) {
    sel
      .transition()
      .ease(d3.easeLinear)
      .duration((bola) => {
        return (bola.next.originT - bola.current.originT) * 1000;
      })
      .attr("cx", (bola) => {
        return bola.next.origin.x;
      })
      .attr("cy", (bola) => {
        return bola.next.origin.y;
      })
      .on("end", function (bola) {
        try {
          ctx2.beginPath();
          ctx2.moveTo(bola.current.origin.x + 400, bola.current.origin.y + 200);
          ctx2.lineTo(bola.next.origin.x + 400, bola.next.origin.y + 200);
          let x = bola.next.origin.x + 400;
          let y = position(bola.next.origin.y + 190);

          ctx2.lineWidth = 3;
          ctx2.strokeStyle = pintaBolaBlanca ? "white" : "black";
          ctx2.stroke();

          ctx2.font = "40px Calibri"; // Use pixels instead of points
          ctx2.fillStyle = "#000000";

          ctx2.fillText(usedLetterconter % 2 == 0 ? "A" : "B", x, y);
          usedLetterconter++;

          bola.current = bola.next;
          bola.next = decideSiguienteBote(bola.current, escena);
          if (nDeIteracion == iteracion2)
            d3.select(this).call(
              siguienteTransicion2(nDeIteracion, escena, pintaBolaBlanca)
            );
        } catch (e) {
          console.error("2");
        }
      });
  };
}

function position(y) {
  if (y <= 0) {
    return y + 50;
  } else {
    return y;
  }
}

var usedLetterconter = 0;

demuestra2(escena, añadeFondoBunimovich);

////

const svg4 = d3
  .select("#contenedor-periodicidad-4")
  .append("svg")
  .style("position", "absolute")
  .attr("width", 800)
  .attr("height", 400);

const stadiumBg4 = svg4.append("g");

d3.select("#contenedor-periodicidad-4").style("position", "relative");

const canvas4 = d3
  .select("#contenedor-periodicidad-4")
  .append("canvas")
  .style("position", "absolute")
  .attr("width", 800)
  .attr("height", 400)
  .style("opacity", 0.7)
  .node();

d3.select("#contenedor-periodicidad-4")
  .append("div")
  .style("position", "relative")
  .style("height", "420px")
  .style("width", "800px");

const ctx4 = canvas4.getContext("2d");

// Esto indica de donde nacen las bolas
const svgScene4 = svg4.append("g").attr("transform", "translate(400, 200)");

let bolas4 = [];

function añadeBolas4(n, centro, direccion, diferencia, velocidad) {
  const angleScale = d3
    .scaleLinear()
    .domain([0, n])
    .range([-diferencia, diferencia]);
  for (let i = 0; i < n; ++i) {
    const rot = Transform.rotate(angleScale(i));
    bolas4.push({
      current: hazBola(centro, rot.apply(direccion).scale(velocidad)),
    });
  }
}

const canvasBola4 = svgScene4.append("g");

document.getElementById("reinicio-4").addEventListener("click", () => {
  demuestraSensibilidad4(escena, añadeFondoBunimovich);
  c = 0;
});

let iteracion4 = 0;
function demuestraSensibilidad4(escena, backgroundFunction) {
  ++iteracion4;
  bolas4 = [];
  ctx4.clearRect(0, 0, 800, 400);
  console.log(backgroundFunction);

  backgroundFunction(stadiumBg4);

  añadeBolas4(
    1,
    new Vec2(100, 100 * Math.sqrt(2)),
    new Vec2(180, 0),
    0.0,
    velocidad
  );

  const todasLasBolas = canvasBola4.selectAll("*");
  todasLasBolas.transition();
  todasLasBolas.remove();

  const bolasVerdes = canvasBola4
    .selectAll("*")
    .data(bolas4)
    .enter()
    .append("circle")
    .each((ball) => {
      ball.next = decideSiguienteBote(ball.current, escena);
    })
    .attr("cx", function (ball) {
      return ball.current.origin.x;
    })
    .attr("cy", function (ball) {
      return ball.current.origin.y;
    })
    .attr("r", 5)
    .attr("fill", "green");
  bolasVerdes.call(siguienteTransicion4(iteracion4, escena, false));
}

let c = 0;

function siguienteTransicion4(nDeIteracion, escena, pintaBolaBlanca) {
  return function (sel) {
    sel
      .transition()
      .ease(d3.easeLinear)
      .duration((bola) => {
        return (bola.next.originT - bola.current.originT) * 1000;
      })
      .attr("cx", (bola) => {
        return bola.next.origin.x;
      })
      .attr("cy", (bola) => {
        return bola.next.origin.y;
      })
      .on("end", function (bola) {
        try {
          ctx4.beginPath();
          ctx4.moveTo(bola.current.origin.x + 400, bola.current.origin.y + 200);
          ctx4.lineTo(bola.next.origin.x + 400, bola.next.origin.y + 200);
          ctx4.lineWidth = 3;
          ctx4.strokeStyle = pintaBolaBlanca ? "white" : "black";
          ctx4.stroke();
          bola.current = bola.next;
          bola.next = decideSiguienteBote(bola.current, escena);
          if (c++ < 10 && nDeIteracion == iteracion4)
            d3.select(this).call(
              siguienteTransicion4(nDeIteracion, escena, pintaBolaBlanca)
            );
        } catch (e) {
          console.error("");
        }
      });
  };
}

demuestraSensibilidad4(escena, añadeFondoBunimovich);

/////
const svg5 = d3
  .select("#contenedor-periodicidad-5")
  .append("svg")
  .style("position", "absolute")
  .attr("width", 800)
  .attr("height", 400);

const stadiumBg5 = svg5.append("g");

d3.select("#contenedor-periodicidad-5").style("position", "relative");

const canvas5 = d3
  .select("#contenedor-periodicidad-5")
  .append("canvas")
  .style("position", "absolute")
  .attr("width", 800)
  .attr("height", 400)
  .style("opacity", 0.7)
  .node();

d3.select("#contenedor-periodicidad-5")
  .append("div")
  .style("position", "relative")
  .style("height", "420px")
  .style("width", "800px");

const ctx5 = canvas5.getContext("2d");

// Esto indica de donde nacen las bolas
const svgScene5 = svg5.append("g").attr("transform", "translate(400, 200)");

let bolas5 = [];

function añadeBolas5(n, centro, direccion, diferencia, velocidad) {
  const angleScale = d3
    .scaleLinear()
    .domain([0, n])
    .range([-diferencia, diferencia]);
  for (let i = 0; i < n; ++i) {
    const rot = Transform.rotate(angleScale(i));
    bolas5.push({
      current: hazBola(centro, rot.apply(direccion).scale(velocidad)),
    });
  }
}

const canvasBola5 = svgScene5.append("g");

document.getElementById("reinicio-5").addEventListener("click", () => {
  demuestraSensibilidad5(escena, añadeFondoBunimovich);
  c5 = 0;
});

let iteracion5 = 0;
function demuestraSensibilidad5(escena, backgroundFunction) {
  ++iteracion5;
  bolas5 = [];
  ctx5.clearRect(0, 0, 800, 400);
  console.log(backgroundFunction);

  backgroundFunction(stadiumBg5);

  añadeBolas5(1, new Vec2(0, 0), new Vec2(-90, -90), 0.0, velocidad);

  añadeBolas5(1, new Vec2(0, 0), new Vec2(90, 90), 0.0, velocidad);

  const todasLasBolas = canvasBola5.selectAll("*");
  todasLasBolas.transition();
  todasLasBolas.remove();

  const bolasVerdes = canvasBola5
    .selectAll("*")
    .data(bolas5)
    .enter()
    .append("circle")
    .each((ball) => {
      ball.next = decideSiguienteBote(ball.current, escena);
    })
    .attr("cx", function (ball) {
      return ball.current.origin.x;
    })
    .attr("cy", function (ball) {
      return ball.current.origin.y;
    })
    .attr("r", 5)
    .attr("fill", "green");
  bolasVerdes.call(siguienteTransicion5(iteracion5, escena, false));
}

let c5 = 0;

function siguienteTransicion5(nDeIteracion, escena, pintaBolaBlanca) {
  return function (sel) {
    sel
      .transition()
      .ease(d3.easeLinear)
      .duration((bola) => {
        return (bola.next.originT - bola.current.originT) * 1000;
      })
      .attr("cx", (bola) => {
        return bola.next.origin.x;
      })
      .attr("cy", (bola) => {
        return bola.next.origin.y;
      })
      .on("end", function (bola) {
        try {
          ctx5.beginPath();
          ctx5.moveTo(bola.current.origin.x + 400, bola.current.origin.y + 200);
          ctx5.lineTo(bola.next.origin.x + 400, bola.next.origin.y + 200);
          ctx5.lineWidth = 3;
          ctx5.strokeStyle = pintaBolaBlanca ? "white" : "black";
          ctx5.stroke();
          bola.current = bola.next;
          bola.next = decideSiguienteBote(bola.current, escena);
          if (c5++ < 20 && nDeIteracion == iteracion5)
            d3.select(this).call(
              siguienteTransicion5(nDeIteracion, escena, pintaBolaBlanca)
            );
        } catch (e) {
          console.error("");
        }
      });
  };
}

demuestraSensibilidad5(escena, añadeFondoBunimovich);

/// Periodicidad, segunda parte: Densidad
const svg3 = d3
  .select("#contenedor-periodicidad-densidad")
  .append("svg")
  .style("position", "absolute")
  .attr("width", 800)
  .attr("height", 400);

const stadiumBg3 = svg3.append("g");

d3.select("#contenedor-periodicidad-densidad").style("position", "relative");

const canvas3 = d3
  .select("#contenedor-periodicidad-densidad")
  .append("canvas")
  .style("position", "absolute")
  .attr("width", 800)
  .attr("height", 400)
  .style("opacity", 0.7)
  .node();

d3.select("#contenedor-periodicidad-densidad")
  .append("div")
  .style("position", "relative")
  .style("height", "420px")
  .style("width", "800px");

const ctx3 = canvas3.getContext("2d");

// Esto indica de donde nacen las bolas
const svgScene3 = svg3.append("g").attr("transform", "translate(400, 200)");

let bolas3 = [];

function añadeBolas3(n, centro, direccion, diferencia, velocidad) {
  const angleScale = d3
    .scaleLinear()
    .domain([0, n])
    .range([-diferencia, diferencia]);
  for (let i = 0; i < n; ++i) {
    const rot = Transform.rotate(angleScale(i + (0.5 * Math.random() - 0.25)));
    bolas3.push({
      current: hazBola(centro, rot.apply(direccion).scale(velocidad)),
    });
  }
}

const canvasBola3 = svgScene3.append("g");

let iteracion3 = 0;
const velocidad3 = 5;
function demuestra3(escena, backgroundFunction) {
  ++iteracion3;
  bolas3 = [];
  ctx3.clearRect(0, 0, 800, 400);
  console.log(backgroundFunction);

  backgroundFunction(stadiumBg3);

  añadeBolas3(
    1,
    new Vec2(0, 0), // center
    new Vec2(0, 90), //direction
    0.00000000001,
    velocidad3
  );

  const todasLasBolas = canvasBola3.selectAll("*");
  todasLasBolas.transition();
  todasLasBolas.remove();

  const bolasVerdes = canvasBola3
    .selectAll("*")
    .data(bolas3)
    .enter()
    .append("circle")
    .each((ball) => {
      ball.next = decideSiguienteBote(ball.current, escena);
    })
    .attr("cx", function (ball) {
      return ball.current.origin.x;
    })
    .attr("cy", function (ball) {
      return ball.current.origin.y;
    })
    .attr("r", 15)
    .attr("fill", "green");
  bolasVerdes.call(siguienteTransicion3(iteracion3, escena, false));

  const bolasAzules = d3.select(bolasVerdes.nodes()[1]).attr("fill", "blue");
  bolasAzules.call(siguienteTransicion3(iteracion3, escena, true));
}

function siguienteTransicion3(nDeIteracion, escena, pintaBolaBlanca) {
  return function (sel) {
    sel
      .transition()
      .ease(d3.easeLinear)
      .duration((bola) => {
        return (bola.next.originT - bola.current.originT) * 1000;
      })
      .attr("cx", (bola) => {
        return bola.next.origin.x;
      })
      .attr("cy", (bola) => {
        return bola.next.origin.y;
      })
      .on("end", function (bola) {
        try {
          ctx3.beginPath();
          ctx3.moveTo(bola.current.origin.x + 400, bola.current.origin.y + 200);
          ctx3.lineTo(bola.next.origin.x + 400, bola.next.origin.y + 200);
          const x = bola.next.origin.x + 400;
          const y = position(bola.next.origin.y + 190);

          ctx3.lineWidth = 3;
          ctx3.strokeStyle = pintaBolaBlanca ? "white" : "black";
          ctx3.stroke();

          if (previousPointPosition != null && usedLetterconter3 >= 2) {
            ctx3.fillStyle = "#f0f";
            ctx3.clearRect(0, 0, 5, 5);
            ctx3.fillText(alphabet[usedLetterconter3 - 2], x, y);
            //console.log("entra aqui");
          }

          ctx3.font = "40px Calibri"; // Use pixels instead of points
          ctx3.fillStyle = "#000000";
          ctx3.fillText(abecedary(), x, y);

          updatepreviousPointPosition(x, y);
          const li = document.createElement("li");
          li.appendChild(document.createTextNode(x, y));
          document.getElementById("puntos-colision").appendChild(li);

          bola.current = bola.next;
          bola.next = decideSiguienteBote(bola.current, escena);
          if (nDeIteracion == iteracion3)
            d3.select(this).call(
              siguienteTransicion3(nDeIteracion, escena, pintaBolaBlanca)
            );
        } catch (e) {
          console.error("3");
        }
      });
  };
}

function updatepreviousPointPosition(x, y) {
  previousPointPosition = [];
  previousPointPosition.push(x);
  previousPointPosition.push(y);
}

const alphabet = "abcdefghijklmnopqrstuvwxyz".toUpperCase().split("");
var usedLetterconter3 = 0;
var alphabetLetter = "A";
var previousPointPosition = null;

function abecedary() {
  //console.log(alphabet.length);
  alphabetLetter = alphabet[usedLetterconter3];
  if (usedLetterconter3 < alphabet.length - 1) {
    usedLetterconter3++;
  } else {
    usedLetterconter3 = 0;
  }

  return alphabetLetter;
}

demuestra3(escena, añadeFondoBunimovich);

/// Transitividad
const svg6 = d3
  .select("#contenedor-transitividad")
  .append("svg")
  .style("position", "absolute")
  .attr("width", 800)
  .attr("height", 400);

const stadiumBg6 = svg6.append("g");

d3.select("#contenedor-transitividad").style("position", "relative");

const canvas6 = d3
  .select("#contenedor-transitividad")
  .append("canvas")
  .style("position", "absolute")
  .attr("width", 800)
  .attr("height", 400)
  .style("opacity", 0.7)
  .node();

d3.select("#contenedor-transitividad")
  .append("div")
  .style("position", "relative")
  .style("height", "420px")
  .style("width", "800px");

const ctx6 = canvas6.getContext("2d");

// Esto indica de donde nacen las bolas
const svgScene6 = svg6.append("g").attr("transform", "translate(400, 200)");

let bolas6 = [];

function añadebolas6(n, centro, direccion, diferencia, velocidad) {
  const angleScale = d3
    .scaleLinear()
    .domain([0, n])
    .range([-diferencia, diferencia]);
  for (let i = 0; i < n; ++i) {
    const rot = Transform.rotate(angleScale(i + (0.5 * Math.random() - 0.25)));
    bolas6.push({
      current: hazBola(centro, rot.apply(direccion).scale(velocidad)),
    });
  }
}

const canvasBola6 = svgScene6.append("g");

let iteracion6 = 0;
function demuestra6(escena, backgroundFunction) {
  contadorIterTransBolaVerde = 0;
  contadorIterTransBolaAzul = 0;
  contadorIterTransBolaAmarillo = 0;
  ++iteracion6;
  bolas6 = [];
  ctx6.clearRect(0, 0, 800, 400);
  console.log(backgroundFunction);

  backgroundFunction(stadiumBg6);

  añadebolas6(
    1,
    new Vec2(147, 190), // center
    new Vec2(30, 52), //direction
    0.0,
    5
  );

  añadebolas6(
    1,
    new Vec2(50, 20), // center
    new Vec2(30, 52.15), //direction
    0.0,
    5
  );

  añadebolas6(
    1,
    new Vec2(-300, -150), // center
    new Vec2(-62, -20), //direction
    0.0,
    5
  );

  const todasLasBolas = canvasBola6.selectAll("*");
  todasLasBolas.transition();
  todasLasBolas.remove();

  const bolasVerdes = canvasBola6
    .selectAll("*")
    .data(bolas6)
    .enter()
    .append("circle")
    .each((ball) => {
      ball.next = decideSiguienteBote(ball.current, escena);
    })
    .attr("cx", function (ball) {
      return ball.current.origin.x;
    })
    .attr("cy", function (ball) {
      return ball.current.origin.y;
    })
    .attr("r", 15)
    .attr("fill", "green");

  const bolasAzules = d3.select(bolasVerdes.nodes()[1]).attr("fill", "blue");
  const bolasAmarillas = d3
    .select(bolasVerdes.nodes()[2])
    .attr("fill", "yellow");

  bolasVerdes.call(
    siguienteTransicion6(
      iteracion6,
      escena,
      true,
      contadorIterTransBolaVerde,
      stopITerBolaVerde,
      "green"
    )
  );

  bolasAzules.call(
    siguienteTransicion6(
      iteracion6,
      escena,
      false,
      contadorIterTransBolaAzul,
      stopITerBolaAzul,
      "blue"
    )
  );

  bolasAmarillas.call(
    siguienteTransicion6(
      iteracion6,
      escena,
      true,
      contadorIterTransBolaAmarillo,
      stopITerBolaAmarillo,
      "yellow"
    )
  );
}

function siguienteTransicion6(
  nDeIteracion,
  escena,
  pintaBolaBlanca,
  contadorIterTransBola,
  stopITerBola,
  color
) {
  //console.log(nDeIteracion);

  return function (sel) {
    sel
      .transition()
      .ease(d3.easeLinear)
      .duration((bola) => {
        return (bola.next.originT - bola.current.originT) * 1000;
      })
      .attr("cx", (bola) => {
        return bola.next.origin.x;
      })
      .attr("cy", (bola) => {
        return bola.next.origin.y;
      })
      .on("end", function (bola) {
        try {
          ctx6.beginPath();
          ctx6.moveTo(bola.current.origin.x + 400, bola.current.origin.y + 200);
          ctx6.lineTo(bola.next.origin.x + 400, bola.next.origin.y + 200);
          let x = bola.next.origin.x + 400;
          let y = bola.next.origin.y + 190;
          console.log(color + ": " + x, y);

          ctx6.lineWidth = 3;
          ctx6.strokeStyle = pintaBolaBlanca ? "white" : "black";
          ctx6.stroke();

          ctx6.font = "18px Calibri"; // Use pixels instead of points
          ctx6.fillStyle = color;
          //console.log(x,y);

          ctx6.strokeStyle = "black";
          ctx6.lineWidth = 2;

          let newPos = getNewPos(
            Math.round((x + Number.EPSILON) * 100) / 100,
            Math.round((y + Number.EPSILON) * 100) / 100,
            color
          );
          console.log("NewPos: " + newPos[0] + ", " + newPos[1]);
          ctx6.strokeText(
            getLetter(color)[contadorIterTransBola] +
              "(" +
              newPos[0].toString() +
              "," +
              newPos[1].toString() +
              ")",
            newPos[0],
            newPos[1]
          );
          ctx6.fillText(
            getLetter(color)[contadorIterTransBola] +
              "(" +
              newPos[0].toString() +
              "," +
              newPos[1].toString() +
              ")",
            newPos[0],
            newPos[1]
          );

          bola.current = bola.next;
          bola.next = decideSiguienteBote(bola.current, escena);

          if (
            nDeIteracion == iteracion6 &&
            contadorIterTransBola++ < stopITerBola
          ) {
            d3.select(this).call(
              siguienteTransicion6(
                nDeIteracion,
                escena,
                pintaBolaBlanca,
                contadorIterTransBola,
                stopITerBola,
                color
              )
            );
            console.log("Contador Iter Transividad " + contadorIterTransBola);
          }
        } catch (e) {
          console.error("2");
        }
      });
  };
}

function getLetter(color) {
  let res = ab;
  switch (color) {
    case "green":
      res = ab;
      break;

    case "blue":
      res = acde;
      break;

    default:
      res = fd;
      break;
  }
  return res;
}

const ab = ["A", "B"];
const acde = ["A", "C", "D", "E"];
const fd = ["F", "D", "C"];

function getNewPos(x, y, color) {
  let res = [x, y];

  if (x > 600) {
    res[0] = x - 100;
  }
  if (color == "blue" && y < 101) {
    res[1] = y + 40;
  } else if (color == "blue" && y > 100) {
    res[1] = y - 20;
  }
  return res;
}

let contadorIterTransBolaVerde = 0;
let contadorIterTransBolaAzul = 0;
let contadorIterTransBolaAmarillo = 0;
const stopITerBolaVerde = 1;
const stopITerBolaAzul = 3;
const stopITerBolaAmarillo = 1;

document
  .getElementById("reinicio-transitividad")
  .addEventListener("click", () => {
    demuestra6(escena, añadeFondoBunimovich);
  });

setTimeout(() => {
  demuestra6(escena, añadeFondoBunimovich);
}, 4000);
