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

function lineLineIntersect(p1, p2, p3, p4) {
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

function makeBall(origin, velocity, originT) {
  return {
    originT: originT || 0,
    origin: origin,
    velocity: velocity,
  };
}

function shapeBounce(shape, ball) {
  const intersection = shape.intersect(ball);
  if (
    intersection.t === Infinity ||
    isNaN(intersection.t) ||
    intersection.t <= ball.originT
  )
    return undefined;

  const normal = shape.normalAt(intersection.point);
  if (normal.dot(ball.velocity.unit()) >= 0) return undefined;

  return makeBall(
    intersection.point,
    ball.velocity.reflect(normal).neg(),
    intersection.t
  );
}

function halfPlane(point, normal) {
  normal = normal.unit();
  const tangent = new Vec2(normal.y, -normal.x);
  const p1 = point.plus(tangent),
    p2 = point.minus(tangent);

  return {
    intersect(ball) {
      const result = lineLineIntersect(
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

function halfInnerCircle(center, radius, normal) {
  const ellipse = new Ellipse()
    .transform(Transform.scale(radius, radius))
    .transform(Transform.translate(center.x, center.y));

  normal = normal.unit();
  return {
    intersect: function (ball) {
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
    normalAt: function (point) {
      return ellipse.gradient.apply(point).unit().neg();
    },
  };
}

//////////////////////////////////////////////////////////////////////////////

function decideNextBounce(ball, scene) {
  const bounces = [];
  scene.forEach(function (object) {
    const intersect = object.intersect(ball);
    if (intersect === undefined) return;
    if (intersect.t === Infinity) return;
    const bounce = shapeBounce(object, ball);
    if (bounce !== undefined) {
      bounces.push(bounce);
    }
  });
  bounces.sort(function (b1, b2) {
    if (b1.originT < b2.originT) return -1;
    else if (b1.originT > b2.originT) return 1;
    else return 0;
  });
  if (bounces.length === 0) {
    throw new Error("ball escaped?");
  }
  return bounces[0];
}

//////////////////////////////////////////////////////////////////////////////

const svg = d3
  .select("#main")
  .append("svg")
  .style("position", "absolute")
  .attr("width", 800)
  .attr("height", 400);

const bunimovichScene = [
  halfPlane(new Vec2(0, 200), new Vec2(0, -1)),
  halfPlane(new Vec2(0, -200), new Vec2(0, 1)),
  halfInnerCircle(new Vec2(-200, 0), 200, new Vec2(1, 0)),
  halfInnerCircle(new Vec2(200, 0), 200, new Vec2(-1, 0)),
];

function addBunimovichBG(sel) {
  sel
    .attr("transform", "translate(400, 200)")
    .append("path")
    .attr(
      "d",
      "M -200 -200 l 400 0 a 200 200 0 0 1 0 400 l -400 0 a 200 200 0 0 1 0 -400 "
    )
    .attr("fill", d3.lab(80, 0, 0));
}

const stadiumBg = svg.append("g");

d3.select("#main").style("position", "relative");

const stadiumTrace = d3
  .select("#main")
  .append("canvas")
  .style("position", "absolute")
  .attr("width", 800)
  .attr("height", 400)
  .style("opacity", 0.7)
  .node();

d3.select("#main")
  .append("div")
  .style("position", "relative")
  .style("height", "420px")
  .style("width", "800px");

const ctx = stadiumTrace.getContext("2d");

// Esto indica de donde nacen las bolas
const svgScene = svg.append("g").attr("transform", "translate(400, 200)");

let balls = [];

function addBalls(n, center, direction, spread, speed) {
  const angleScale = d3.scaleLinear().domain([0, n]).range([-spread, spread]);
  for (let i = 0; i < n; ++i) {
    const rot = Transform.rotate(angleScale(i + (0.5 * Math.random() - 0.25)));
    balls.push({
      current: makeBall(center, rot.apply(direction).scale(speed)),
    });
  }
}

const ballPen = svgScene.append("g");

const [selector, scene, backgroundFunction] = [
  "#reset",
  bunimovichScene,
  addBunimovichBG,
];

d3.select(selector)
  .append("span")
  .style("min-width", "200px")
  .style("display", "inline-block")
  .text("Scene: " + selector.substring(7));

d3.select(selector)
  .append("button")
  .text("300 balls")
  .on("click", function () {
    reset(300, scene, backgroundFunction);
  });

d3.select(selector).append("span").text(" ");

d3.select(selector)
  .append("button")
  .text("1 ball")
  .on("click", function () {
    reset(1, scene, backgroundFunction);
  });

d3.select("#clear-tracks")
  .append("button")
  .text("clear tracks")
  .on("click", function () {
    ctx.clearRect(0, 0, 800, 400);
  });

let batch = 0;
const speed = 500;
function reset(howMany, scene, backgroundFunction) {
  ++batch;
  balls = [];
  ctx.clearRect(0, 0, 800, 400);

  backgroundFunction(stadiumBg);

  const randomAngle = Math.random() * 2 * Math.PI;

  addBalls(
    howMany,
    new Vec2(Math.random() * 100 - 50, Math.random() * 100 - 50),
    new Vec2(Math.cos(randomAngle), Math.sin(randomAngle)),
    Math.random() * 0.25,
    speed
  );

  const all = ballPen.selectAll("*");
  all.transition();
  all.remove();

  ballPen
    .selectAll("*")
    .data(balls)
    .enter()
    .append("circle")
    .each(function (ball) {
      ball.next = decideNextBounce(ball.current, scene);
    })
    .attr("cx", function (ball) {
      return ball.current.origin.x;
    })
    .attr("cy", function (ball) {
      return ball.current.origin.y;
    })
    .attr("r", 2)
    .attr("fill", "peru")
    .call(transitionToNext(batch, scene));
}

reset(300, bunimovichScene, addBunimovichBG);

function transitionToNext(whichBatch, scene) {
  return function (sel) {
    sel
      .transition()
      .ease(d3.easeLinear)
      .duration(function (ball) {
        return (ball.next.originT - ball.current.originT) * 1000;
      })
      .attr("cx", function (ball) {
        return ball.next.origin.x;
      })
      .attr("cy", function (ball) {
        return ball.next.origin.y;
      })
      .on("end", function (ball) {
        ctx.beginPath();
        ctx.moveTo(ball.current.origin.x + 400, ball.current.origin.y + 200);
        ctx.lineTo(ball.next.origin.x + 400, ball.next.origin.y + 200);
        ctx.strokeStyle = "blue";
        ctx.stroke();
        ball.current = ball.next;
        ball.next = decideNextBounce(ball.current, scene);
        if (whichBatch === batch) {
          d3.select(this).call(transitionToNext(whichBatch, scene));
        }
      });
  };
}
