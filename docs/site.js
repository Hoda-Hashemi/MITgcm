document.body.classList.add("js-ready");

const reducedMotionQuery = window.matchMedia("(prefers-reduced-motion: reduce)");
const prefersReducedMotion = reducedMotionQuery.matches;

const navToggle = document.querySelector("[data-nav-toggle]");
const navDrawer = document.querySelector("[data-nav-drawer]");
const navClose = document.querySelector("[data-nav-close]");
const navBackdrop = document.querySelector("[data-nav-backdrop]");
const navLinks = Array.from(document.querySelectorAll(".side-nav a[href^='#']"));
const sections = navLinks
  .map((link) => document.querySelector(link.getAttribute("href")))
  .filter(Boolean);

function setNavOpen(open) {
  document.body.classList.toggle("nav-open", open);
  if (navToggle) {
    navToggle.setAttribute("aria-expanded", String(open));
  }
}

function setActiveNav(hash) {
  let activeLink = null;
  navLinks.forEach((link) => {
    const active = link.getAttribute("href") === hash;
    link.classList.toggle("is-active", active);
    link.classList.toggle("active", active);
    if (active) {
      activeLink = link;
    }
  });

  document.querySelectorAll(".nav-group").forEach((group) => {
    const hasActive = activeLink ? group.contains(activeLink) : false;
    group.classList.toggle("has-active", hasActive);
    if (hasActive) {
      group.open = true;
    }
  });
}

function scrollToSection(hash) {
  const target = document.querySelector(hash);
  if (!target) {
    return;
  }
  target.scrollIntoView({
    behavior: prefersReducedMotion ? "auto" : "smooth",
    block: "start",
  });
  window.history.replaceState(null, "", hash);
  setActiveNav(hash);
}

if (navToggle) {
  navToggle.addEventListener("click", () => setNavOpen(!document.body.classList.contains("nav-open")));
}

[navClose, navBackdrop].forEach((control) => {
  if (control) {
    control.addEventListener("click", () => setNavOpen(false));
  }
});

navLinks.forEach((link) => {
  link.addEventListener("click", (event) => {
    const hash = link.getAttribute("href");
    if (!hash || hash === "#") {
      return;
    }
    event.preventDefault();
    scrollToSection(hash);
    setNavOpen(false);
  });
});

let activePending = false;

function updateActiveFromScroll() {
  const anchorOffset = Math.min(180, window.innerHeight * 0.28);
  let current = sections[0];
  sections.forEach((section) => {
    if (section.getBoundingClientRect().top <= anchorOffset) {
      current = section;
    }
  });
  if (current) {
    setActiveNav(`#${current.id}`);
  }
  activePending = false;
}

function requestActiveUpdate() {
  if (activePending) {
    return;
  }
  activePending = true;
  window.requestAnimationFrame(updateActiveFromScroll);
}

if (sections.length) {
  updateActiveFromScroll();
  window.addEventListener("scroll", requestActiveUpdate, { passive: true });
  window.addEventListener("resize", requestActiveUpdate);
}

const revealTargets = document.querySelectorAll(
  ".section-card, .description-block, .metric-card, .subfigure"
);

if ("IntersectionObserver" in window) {
  const revealObserver = new IntersectionObserver(
    (entries, observer) => {
      entries.forEach((entry) => {
        if (!entry.isIntersecting) {
          return;
        }
        entry.target.classList.add("is-visible");
        observer.unobserve(entry.target);
      });
    },
    { rootMargin: "0px 0px -8% 0px", threshold: 0.08 }
  );

  revealTargets.forEach((target) => revealObserver.observe(target));
} else {
  revealTargets.forEach((target) => target.classList.add("is-visible"));
}

const parallaxLayers = Array.from(document.querySelectorAll("[data-parallax]"));
let parallaxPending = false;

function updateParallax() {
  const scrollY = window.scrollY || window.pageYOffset || 0;
  parallaxLayers.forEach((layer) => {
    const speed = Number(layer.getAttribute("data-parallax")) || 0;
    layer.style.transform = `translate3d(0, ${scrollY * speed}px, 0)`;
  });
  parallaxPending = false;
}

function requestParallaxUpdate() {
  if (parallaxPending) {
    return;
  }
  parallaxPending = true;
  window.requestAnimationFrame(updateParallax);
}

if (parallaxLayers.length && !prefersReducedMotion) {
  updateParallax();
  window.addEventListener("scroll", requestParallaxUpdate, { passive: true });
}

function closeModal(modal) {
  const image = modal.querySelector("[data-modal-img]");
  modal.hidden = true;
  if (image) {
    image.removeAttribute("src");
  }
}

document.addEventListener("click", (event) => {
  const tab = event.target.closest("[data-tab-target]");
  if (tab) {
    const tabList = tab.closest(".tab-list");
    const tabs = tabList ? tabList.querySelectorAll("[data-tab-target]") : [];
    const tabsRoot = tab.closest(".field-tabs");
    const panels = tabsRoot ? tabsRoot.querySelectorAll("[data-tab-panel]") : [];
    const target = tab.getAttribute("data-tab-target");

    tabs.forEach((button) => {
      const active = button === tab;
      button.classList.toggle("is-active", active);
      button.setAttribute("aria-selected", String(active));
    });

    panels.forEach((panel) => {
      const active = panel.getAttribute("data-tab-panel") === target;
      panel.classList.toggle("is-active", active);
      panel.hidden = !active;
    });

    if (tabsRoot && window.MathJax && window.MathJax.typesetPromise) {
      window.MathJax.typesetPromise([tabsRoot]);
    }
    return;
  }

  const openButton = event.target.closest("[data-modal-src]");
  const modal = document.querySelector("[data-modal]");
  if (openButton && modal) {
    const image = modal.querySelector("[data-modal-img]");
    const caption = modal.querySelector("[data-modal-caption]");
    image.src = openButton.getAttribute("data-modal-src");
    image.alt = openButton.getAttribute("data-modal-caption") || "";
    caption.textContent = openButton.getAttribute("data-modal-caption") || "";
    modal.hidden = false;
    return;
  }

  if (modal && (event.target.closest("[data-modal-close]") || event.target === modal)) {
    closeModal(modal);
  }
});

document.addEventListener("keydown", (event) => {
  if (event.key !== "Escape") {
    return;
  }

  if (document.body.classList.contains("nav-open")) {
    setNavOpen(false);
    return;
  }

  const modal = document.querySelector("[data-modal]");
  if (modal && !modal.hidden) {
    closeModal(modal);
  }
});
